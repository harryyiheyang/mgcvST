#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include "liu_tail.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// PCAlearning approximate Liu path.
//
// Symmetric q x q matrices are stored in weighted-vech form: the upper
// triangle in column-major order (entry (i, j), i <= j, at j (j + 1) / 2 + i)
// with off-diagonal weight sqrt(2), so that inner products of weighted-vech
// vectors are Frobenius inner products. Training matrices H_j are kept in this
// form as float32 raw vectors; the basis B is returned in this form (L x r,
// L = q (q + 1) / 2) and its columns B_k are the symmetric basis matrices.
//
// Trace tables: with X = sum_k x_k B_k and Y = sum_k y_k B_k,
// tr((XY)^s) = sum kappa_s(x)' Tsym_s kappa_s(y), where kappa_s are the degree-s
// monomials indexed by nondecreasing index tuples in lexicographic order.
//   s = 1: Tsym1 = B'B.
//   s = 2: tr(XYXY) = sum tr(Y_ab Y_cd), Y_ab = B_a B_b, and tr(Y_ab Y_cd) =
//          <Y_ba, Y_cd>: one Gram matrix of the r^2 products.
//   s = 3: tr((XY)^3) = <XYX, YXY>. XYX has symmetric coefficients
//          N_{m,c} = sum_{perm m} B_m1 B_c B_m2 (m a 2-multiset); Tsym3 is a
//          scatter of the Gram matrix of the N_{m,c}.
//   s = 4: tr((XY)^4) = tr(E^2), E = XYXY with coefficients W_{m,n}
//          (m, n 2-multisets) satisfying W_{m,n}' = W_{n,m}, so that
//          tr(W_u W_v) = <W_swap(u), W_v>. The W_{m,n} are formed on a y grid
//          (y = e_b and e_b + e_c) as N^(y)_m Y and one GEMM per grid point, in
//          row blocks restricted to the upper triangle; Tsym4 is a scatter of
//          the resulting Gram matrix.
// The GEMMs of levels 2-4 run in float32 with double accumulation across row
// blocks, as in the validated experiment.

namespace {

using MatD = Eigen::MatrixXd;
using MatF = Eigen::MatrixXf;
using Index = Eigen::Index;
using mgcvst_liu::liu_log_p;

double now() {
  return std::chrono::duration<double>(
    std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Multisets of size s over {0, ..., r - 1}: nondecreasing tuples in
// lexicographic order (first coordinate slowest).
struct Multisets {
  int r = 0, s = 0, d = 0;
  std::vector<int> tuple;  // d x s, row-major
  std::vector<int> id;     // lexicographic index of a sorted tuple -> multiset id
  int at(int t, int k) const { return tuple[(size_t)t * s + k]; }
  int find(int* x) const {
    std::sort(x, x + s);
    long z = 0;
    for (int k = 0; k < s; ++k) z = z * r + x[k];
    return id[z];
  }
};

Multisets make_multisets(int r, int s) {
  Multisets m;
  m.r = r;
  m.s = s;
  long N = 1;
  for (int k = 0; k < s; ++k) N *= r;
  m.id.assign(N, -1);
  std::vector<int> x(s);
  for (long t = 0; t < N; ++t) {
    long z = t;
    for (int k = s - 1; k >= 0; --k) {
      x[k] = z % r;
      z /= r;
    }
    bool ok = true;
    for (int k = 1; k < s; ++k) if (x[k] < x[k - 1]) ok = false;
    if (!ok) continue;
    m.id[t] = m.d++;
    for (int k = 0; k < s; ++k) m.tuple.push_back(x[k]);
  }
  return m;
}

Rcpp::IntegerMatrix multiset_matrix(const Multisets& m) {
  Rcpp::IntegerMatrix out(m.d, m.s);
  for (int t = 0; t < m.d; ++t) for (int k = 0; k < m.s; ++k) out(t, k) = m.at(t, k) + 1;
  return out;
}

// Monomials kappa_s(c_g) for genes g in cols (d x |cols|).
MatD monomials(const Multisets& m, const MatD& C, const std::vector<int>& cols) {
  MatD K(m.d, cols.size());
  for (size_t g = 0; g < cols.size(); ++g) {
    for (int t = 0; t < m.d; ++t) {
      double z = 1;
      for (int k = 0; k < m.s; ++k) z *= C(cols[g], m.at(t, k));
      K(t, g) = z;
    }
  }
  return K;
}

// G(upper tiles) += P' P for P (len x n) float; tiles run in parallel with
// float products and double accumulation.
void gram_accumulate(const float* data, Index len, Index n, MatD& G,
                     int threads, int tile) {
  Eigen::Map<const MatF> P(data, len, n);
  std::vector<std::pair<Index, Index> > tasks;
  for (Index u = 0; u < n; u += tile) for (Index v = u; v < n; v += tile) tasks.push_back({u, v});
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (long k = 0; k < (long)tasks.size(); ++k) {
    const Index u = tasks[k].first, v = tasks[k].second;
    const Index nu = std::min<Index>(tile, n - u), nv = std::min<Index>(tile, n - v);
    MatF tmp = P.middleCols(u, nu).transpose() * P.middleCols(v, nv);
    G.block(u, v, nu, nv) += tmp.cast<double>();
  }
}

void symmetrize_upper(MatD& G) {
  for (Index j = 0; j < G.cols(); ++j) for (Index i = j + 1; i < G.rows(); ++i) G(i, j) = G(j, i);
}

// Symmetric q x q matrix from a full (q^2) or weighted-vech (q (q + 1) / 2) column.
MatD unpack_symmetric(const double* x, int q, bool full) {
  MatD M(q, q);
  if (full) {
    Eigen::Map<const MatD> F(x, q, q);
    M = 0.5 * (F + F.transpose());
    return M;
  }
  const double w = 1 / std::sqrt(2.0);
  Index p = 0;
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < j; ++i) M(i, j) = M(j, i) = w * x[p++];
    M(j, j) = x[p++];
  }
  return M;
}

template <class T>
void pack_symmetric(const MatD& M, T* out) {
  const double w = std::sqrt(2.0);
  const int q = M.rows();
  Index p = 0;
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < j; ++i) out[p++] = (T)(w * 0.5 * (M(i, j) + M(j, i)));
    out[p++] = (T)M(j, j);
  }
}

struct Tables {
  int r = 0;
  Multisets ms[5];
  MatD T[5];
};

Tables parse_tables(const Rcpp::List& tables) {
  Tables t;
  Rcpp::NumericMatrix T1 = tables["Tsym1"];
  t.r = T1.nrow();
  for (int s = 1; s <= 4; ++s) {
    t.ms[s] = make_multisets(t.r, s);
    Rcpp::NumericMatrix Ts = tables["Tsym" + std::to_string(s)];
    t.T[s] = Eigen::Map<MatD>(Ts.begin(), Ts.nrow(), Ts.ncol());
  }
  return t;
}

const float* raw_float(SEXP x) { return reinterpret_cast<const float*>(RAW(x)); }

Rcpp::CharacterVector pair_columns(bool moments) {
  if (moments) return Rcpp::CharacterVector::create("U", "t1", "t2", "t3", "t4",
    "logp_two_sided", "logp_positive", "logp_negative");
  return Rcpp::CharacterVector::create("U", "logp_two_sided", "logp_positive", "logp_negative");
}

inline void write_pair(double* out, Index rows, Index row, bool moments, double U,
                       double t1, double t2, double t3, double t4) {
  double lp[3];
  liu_log_p(U, t1, t2, t3, t4, lp);
  out[row] = U;
  int k = 1;
  if (moments) {
    out[row + rows] = t1;
    out[row + 2 * rows] = t2;
    out[row + 3 * rows] = t3;
    out[row + 4 * rows] = t4;
    k = 5;
  }
  for (int z = 0; z < 3; ++z) out[row + (k + z) * rows] = lp[z];
}

} // namespace

// Natural-log Liu p-values (two-sided, positive, negative) of signed scores U
// from trace moments t1..t4 = tr((H_i H_j)^s).
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_liu_logp_cpp(const Rcpp::NumericVector& U,
                                        const Rcpp::NumericVector& t1,
                                        const Rcpp::NumericVector& t2,
                                        const Rcpp::NumericVector& t3,
                                        const Rcpp::NumericVector& t4,
                                        int threads = 1) {
  const R_xlen_t n = U.size();
  Rcpp::NumericMatrix out(n, 3);
  double* o = out.begin();
  const double *u = U.begin(), *a = t1.begin(), *b = t2.begin(), *c = t3.begin(), *d = t4.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1024)
#endif
  for (R_xlen_t k = 0; k < n; ++k) {
    double lp[3];
    liu_log_p(u[k], a[k], b[k], c[k], d[k], lp);
    o[k] = lp[0];
    o[k + n] = lp[1];
    o[k + 2 * n] = lp[2];
  }
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("two_sided", "positive", "negative");
  return out;
}

// Weighted-vech float32 packing of a symmetric matrix (symmetrized first).
// [[Rcpp::export]]
Rcpp::RawVector mgcvst_pca_pack_cpp(const Eigen::Map<Eigen::MatrixXd> M) {
  const Index q = M.rows(), L = q * (q + 1) / 2;
  Rcpp::RawVector out(4 * L);
  pack_symmetric<float>(MatD(M), reinterpret_cast<float*>(RAW(out)));
  return out;
}

// Gram matrix G_jk = tau_j tau_k <H_j, H_k>_F of packed float32 training
// matrices, accumulated in double over row chunks.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_gram_cpp(const Rcpp::List& packed,
                                        const Rcpp::NumericVector& tau,
                                        int threads = 1, int chunk = 8192) {
  const int n = packed.size();
  const Index L = XLENGTH(packed[0]) / 4;
  std::vector<const float*> ptr(n);
  for (int j = 0; j < n; ++j) ptr[j] = raw_float(packed[j]);
  const Index nb = (L + chunk - 1) / chunk;
  std::vector<MatD> acc(threads, MatD::Zero(n, n));
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
#ifdef _OPENMP
    const int t = omp_get_thread_num();
#else
    const int t = 0;
#endif
    MatD V(chunk, n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (Index b = 0; b < nb; ++b) {
      const Index r0 = b * chunk, m = std::min<Index>(chunk, L - r0);
      for (int j = 0; j < n; ++j) {
        const float* x = ptr[j] + r0;
        for (Index i = 0; i < m; ++i) V(i, j) = x[i];
      }
      acc[t].selfadjointView<Eigen::Lower>().rankUpdate(V.topRows(m).transpose());
    }
  }
  MatD G = MatD::Zero(n, n);
  for (auto& a : acc) G += a;
  for (int j = 0; j < n; ++j) for (int i = j; i < n; ++i) G(j, i) = G(i, j) = G(i, j) * tau[i] * tau[j];
  return Rcpp::wrap(G);
}

// Basis in weighted-vech form: B = V R, V = [tau_j vech_w(H_j)], R = n x r
// (eigenvectors scaled by 1 / sqrt(eigenvalues)).
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_basis_cpp(const Rcpp::List& packed,
                                         const Rcpp::NumericVector& tau,
                                         const Eigen::Map<Eigen::MatrixXd> R,
                                         int threads = 1, int chunk = 8192) {
  const int n = packed.size(), r = R.cols();
  const Index L = XLENGTH(packed[0]) / 4;
  std::vector<const float*> ptr(n);
  for (int j = 0; j < n; ++j) ptr[j] = raw_float(packed[j]);
  MatD TR = Eigen::Map<const Eigen::VectorXd>(tau.begin(), n).asDiagonal() * R;
  Rcpp::NumericMatrix out(L, r);
  Eigen::Map<MatD> B(out.begin(), L, r);
  const Index nb = (L + chunk - 1) / chunk;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    MatD V(chunk, n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (Index b = 0; b < nb; ++b) {
      const Index r0 = b * chunk, m = std::min<Index>(chunk, L - r0);
      for (int j = 0; j < n; ++j) {
        const float* x = ptr[j] + r0;
        for (Index i = 0; i < m; ++i) V(i, j) = x[i];
      }
      B.middleRows(r0, m).noalias() = V.topRows(m) * TR;
    }
  }
  return out;
}

// Trace tables Tsym1..Tsym4 of a symmetric basis B given as q^2 x r (full
// columns, symmetrized) or q (q + 1) / 2 x r (weighted vech).
// [[Rcpp::export]]
Rcpp::List mgcvst_pca_tables_cpp(const Rcpp::NumericMatrix& B, int q,
                                 int threads = 1, int block = 32, int tile = 192) {
  const double t_start = now();
  const int r = B.ncol();
  const Index qq = (Index)q * q;
  const bool full = B.nrow() == qq;
  const Multisets m1 = make_multisets(r, 1), m2 = make_multisets(r, 2),
    m3 = make_multisets(r, 3), m4 = make_multisets(r, 4);
  const int d2 = m2.d;
  std::vector<double> timing;
  std::vector<std::string> tnames;
  std::vector<double> memory;
  std::vector<std::string> mnames;
  auto mark = [&](const char* name, double t0) { timing.push_back(now() - t0); tnames.push_back(name); };

  // B_k (double and float) and Y_ab = B_a B_b (double GEMM, stored float).
  double t0 = now();
  std::vector<MatD> Bd(r);
  std::vector<MatF> Bf(r);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int k = 0; k < r; ++k) {
    Bd[k] = unpack_symmetric(&B(0, k), q, full);
    Bf[k] = Bd[k].cast<float>();
  }
  std::vector<MatF> Yf(r * r);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int t = 0; t < d2; ++t) {
    const int a = m2.at(t, 0), b = m2.at(t, 1);
    MatD y = Bd[a] * Bd[b];
    Yf[a * r + b] = y.cast<float>();
    if (a != b) Yf[b * r + a] = y.transpose().cast<float>();
  }
  MatD T1(r, r);
  for (int a = 0; a < r; ++a) for (int b = 0; b < r; ++b) T1(a, b) = (Bd[a].array() * Bd[b].array()).sum();
  Bd.clear();
  mark("products", t0);
  memory.push_back(4.0 * qq * (r + r * r)); mnames.push_back("B_Y");

  // Level 2: Gram of the r^2 products over column blocks of full entries.
  t0 = now();
  MatD Gy = MatD::Zero(r * r, r * r);
  {
    MatF P((Index)q * block, r * r);
    for (int c0 = 0; c0 < q; c0 += block) {
      const int w = std::min(block, q - c0);
      const Index len = (Index)q * w;
      for (int t = 0; t < r * r; ++t) std::memcpy(P.data() + (Index)t * len, Yf[t].data() + (Index)c0 * q, sizeof(float) * len);
      gram_accumulate(P.data(), len, r * r, Gy, threads, tile);
    }
  }
  symmetrize_upper(Gy);
  MatD T2 = MatD::Zero(d2, d2);
  for (int a1 = 0; a1 < r; ++a1) for (int b1 = 0; b1 < r; ++b1)
    for (int a2 = 0; a2 < r; ++a2) for (int b2 = 0; b2 < r; ++b2) {
      int x[2] = {a1, a2}, y[2] = {b1, b2};
      T2(m2.find(x), m2.find(y)) += Gy(b1 * r + a1, a2 * r + b2);
    }
  mark("level2", t0);

  // Level 3: N_{m,c} (symmetric, float), Gram over weighted upper triangles.
  t0 = now();
  const int nN = d2 * r;
  std::vector<MatF> Nf(nN);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int u = 0; u < nN; ++u) {
    const int m = u / r, c = u % r, a1 = m2.at(m, 0), a2 = m2.at(m, 1);
    MatF P = Bf[a1] * Yf[c * r + a2];
    if (a1 == a2) Nf[u] = 0.5f * (P + P.transpose());
    else Nf[u] = P + P.transpose();
  }
  Yf.clear();
  mark("form3", t0);
  memory.push_back(4.0 * qq * nN); mnames.push_back("N");
  t0 = now();
  MatD G3 = MatD::Zero(nN, nN);
  {
    const float sq2 = (float)std::sqrt(2.0);
    MatF P((Index)q * block, nN);
    for (int c0 = 0; c0 < q; c0 += block) {
      const int w = std::min(block, q - c0), h = c0 + w;
      const Index len = (Index)h * w;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
      for (int u = 0; u < nN; ++u) {
        float* dst = P.data() + (Index)u * len;
        for (int jj = 0; jj < w; ++jj) {
          const int col = c0 + jj;
          const float* src = Nf[u].data() + (Index)col * q;
          float* o = dst + (Index)jj * h;
          for (int i = 0; i < h; ++i) o[i] = i < col ? sq2 * src[i] : (i == col ? src[i] : 0.0f);
        }
      }
      gram_accumulate(P.data(), len, nN, G3, threads, tile);
    }
  }
  symmetrize_upper(G3);
  MatD T3 = MatD::Zero(m3.d, m3.d);
  for (int m = 0; m < d2; ++m) for (int b = 0; b < r; ++b)
    for (int n = 0; n < d2; ++n) for (int a = 0; a < r; ++a) {
      int x[3] = {m2.at(m, 0), m2.at(m, 1), a}, y[3] = {m2.at(n, 0), m2.at(n, 1), b};
      T3(m3.find(x), m3.find(y)) += G3(m * r + b, n * r + a);
    }
  mark("gram3", t0);

  // Level 4: W_{m,n} on the grid y = e_b, e_b + e_c; column u = n d2 + m.
  t0 = now();
  std::vector<MatF> Yg(d2);
  for (int n = 0; n < d2; ++n) {
    const int b = m2.at(n, 0), c = m2.at(n, 1);
    Yg[n] = b == c ? Bf[b] : MatF(Bf[b] + Bf[c]);
  }
  const int nW = d2 * d2;
  MatD Gu = MatD::Zero(nW, nW);
  double t_form4 = 0, t_gram4 = 0;
  {
    MatF P((Index)block * q, nW);
    memory.push_back(4.0 * P.size()); mnames.push_back("W_panel");
    const float isq2 = (float)(1 / std::sqrt(2.0));
    for (int r0 = 0; r0 < q; r0 += block) {
      const int w = std::min(block, q - r0), cw = q - r0;
      const Index len = (Index)w * cw;
      double tf = now();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
      {
        MatF Nbuf(q, (Index)d2 * w);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int n = 0; n < d2; ++n) {
          const int b = m2.at(n, 0), c = m2.at(n, 1);
          for (int m = 0; m < d2; ++m) {
            if (b == c) Nbuf.middleCols((Index)m * w, w) = Nf[m * r + b].middleCols(r0, w);
            else Nbuf.middleCols((Index)m * w, w) = Nf[m * r + b].middleCols(r0, w) + Nf[m * r + c].middleCols(r0, w);
          }
          Eigen::Map<MatF> dest(P.data() + (Index)n * d2 * len, cw, (Index)d2 * w);
          dest.noalias() = Yg[n].middleCols(r0, cw).transpose() * Nbuf;
        }
      }
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
      for (int u = 0; u < nW; ++u) {
        const int n = u / d2, m = u % d2, b = m2.at(n, 0), c = m2.at(n, 1);
        float* x = P.data() + (Index)u * len;
        if (b != c) {
          int bb[2] = {b, b}, cc[2] = {c, c};
          const float* y1 = P.data() + ((Index)m2.find(bb) * d2 + m) * len;
          const float* y2 = P.data() + ((Index)m2.find(cc) * d2 + m) * len;
          for (Index k = 0; k < len; ++k) x[k] -= y1[k] + y2[k];
        }
      }
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
      for (int u = 0; u < nW; ++u) {
        float* x = P.data() + (Index)u * len;
        for (int i = 0; i < w; ++i) {
          float* row = x + (Index)i * cw;
          for (int c = 0; c < i; ++c) row[c] = 0.0f;
          row[i] *= isq2;
        }
      }
      double tg = now();
      t_form4 += tg - tf;
      gram_accumulate(P.data(), len, nW, Gu, threads, tile);
      t_gram4 += now() - tg;
    }
  }
  Nf.clear();
  symmetrize_upper(Gu);
  MatD T4 = MatD::Zero(m4.d, m4.d);
  std::vector<int> sum4((size_t)d2 * d2);
  for (int a = 0; a < d2; ++a) for (int b = 0; b < d2; ++b) {
    int x[4] = {m2.at(a, 0), m2.at(a, 1), m2.at(b, 0), m2.at(b, 1)};
    sum4[(size_t)a * d2 + b] = m4.find(x);
  }
  for (int u = 0; u < nW; ++u) {
    const int n1 = u / d2, k1 = u % d2, su = k1 * d2 + n1;
    for (int v = 0; v < nW; ++v) {
      const int n2 = v / d2, k2 = v % d2, sv = k2 * d2 + n2;
      T4(sum4[(size_t)k1 * d2 + k2], sum4[(size_t)n1 * d2 + n2]) += Gu(su, v) + Gu(u, sv);
    }
  }
  timing.push_back(t_form4); tnames.push_back("form4");
  timing.push_back(t_gram4); tnames.push_back("gram4");
  mark("level4", t0);
  mark("total", t_start);
  memory.push_back(8.0 * nW * nW); mnames.push_back("G4");

  Rcpp::NumericVector tv(timing.begin(), timing.end()), mv(memory.begin(), memory.end());
  tv.names() = tnames;
  mv.names() = mnames;
  return Rcpp::List::create(
    Rcpp::Named("Tsym1") = T1, Rcpp::Named("Tsym2") = T2,
    Rcpp::Named("Tsym3") = T3, Rcpp::Named("Tsym4") = T4,
    Rcpp::Named("ms1") = multiset_matrix(m1), Rcpp::Named("ms2") = multiset_matrix(m2),
    Rcpp::Named("ms3") = multiset_matrix(m3), Rcpp::Named("ms4") = multiset_matrix(m4),
    Rcpp::Named("r") = r, Rcpp::Named("q") = q,
    Rcpp::Named("timing") = tv, Rcpp::Named("memory") = mv);
}

// Approximate pair traces and Liu log p for listed pairs (1-based i, j):
// U = a_i' a_j, t_s = kappa_s(c_i)' Tsym_s kappa_s(c_j) (s = 1: c_i' Tsym1 c_j).
// Columns U, [t1..t4,] logp_two_sided, logp_positive, logp_negative.
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_pairs_cpp(const Eigen::Map<Eigen::MatrixXd> A,
                                         const Eigen::Map<Eigen::MatrixXd> C,
                                         const Rcpp::List& tables,
                                         const Rcpp::IntegerVector& i,
                                         const Rcpp::IntegerVector& j,
                                         int threads = 1, bool moments = true) {
  const Tables tab = parse_tables(tables);
  const int n = C.rows();
  const R_xlen_t np = i.size();
  const MatD Cd = C, CT = Cd * tab.T[1];
  std::vector<int> all(n);
  for (int g = 0; g < n; ++g) all[g] = g;
  std::vector<char> used(n, 0);
  for (R_xlen_t k = 0; k < np; ++k) used[i[k] - 1] = 1;
  std::vector<int> ucols;
  std::vector<int> slot(n, -1);
  for (int g = 0; g < n; ++g) if (used[g]) { slot[g] = ucols.size(); ucols.push_back(g); }
  MatD K[5], S[5];
  for (int s = 2; s <= 4; ++s) {
    K[s] = monomials(tab.ms[s], Cd, all);
    MatD Ku = monomials(tab.ms[s], Cd, ucols);
    S[s].resize(tab.ms[s].d, ucols.size());
    const Index nu = ucols.size(), cb = 64;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (Index b = 0; b < (nu + cb - 1) / cb; ++b) {
      const Index c0 = b * cb, w = std::min(cb, nu - c0);
      S[s].middleCols(c0, w).noalias() = tab.T[s] * Ku.middleCols(c0, w);
    }
  }
  const int ncol = moments ? 8 : 4;
  Rcpp::NumericMatrix out(np, ncol);
  double* o = out.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 256)
#endif
  for (R_xlen_t k = 0; k < np; ++k) {
    const int a = i[k] - 1, b = j[k] - 1, sa = slot[a];
    const double U = A.col(a).dot(A.col(b));
    const double t1 = CT.row(a).dot(Cd.row(b));
    const double t2 = S[2].col(sa).dot(K[2].col(b));
    const double t3 = S[3].col(sa).dot(K[3].col(b));
    const double t4 = S[4].col(sa).dot(K[4].col(b));
    write_pair(o, np, k, moments, U, t1, t2, t3, t4);
  }
  Rcpp::colnames(out) = pair_columns(moments);
  return out;
}

// All pairs (i, j), first <= i <= last, i < j <= n (1-based), ordered by i
// then j. Per block of genes i, traces are GEMMs K_s(i)' Tsym_s K_s(j).
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pca_pairs_block_cpp(const Eigen::Map<Eigen::MatrixXd> A,
                                               const Eigen::Map<Eigen::MatrixXd> C,
                                               const Rcpp::List& tables,
                                               int first, int last,
                                               int threads = 1, bool moments = true,
                                               int gene_block = 32, int pair_block = 1024) {
  const Tables tab = parse_tables(tables);
  const int n = C.rows();
  const int i0 = first - 1, ni = last - first + 1;
  std::vector<Index> offset(ni + 1, 0);
  for (int k = 0; k < ni; ++k) offset[k + 1] = offset[k] + (n - 1 - (i0 + k));
  const Index np = offset[ni];
  const MatD Cd = C, CT = Cd * tab.T[1];
  std::vector<int> all(n);
  for (int g = 0; g < n; ++g) all[g] = g;
  MatD K[5], S[5];
  for (int s = 2; s <= 4; ++s) {
    K[s] = monomials(tab.ms[s], Cd, all);
    MatD Ki = K[s].middleCols(i0, ni);
    S[s].resize(tab.ms[s].d, ni);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (int b = 0; b < (ni + 63) / 64; ++b) {
      const int c0 = b * 64, w = std::min(64, ni - c0);
      S[s].middleCols(c0, w).noalias() = tab.T[s] * Ki.middleCols(c0, w);
    }
  }
  struct Task { int ib, nb, j0, nj; };
  std::vector<Task> tasks;
  for (int ib = 0; ib < ni; ib += gene_block) {
    const int nb = std::min(gene_block, ni - ib);
    for (int j0 = i0 + ib + 1; j0 < n; j0 += pair_block) tasks.push_back({ib, nb, j0, std::min(pair_block, n - j0)});
  }
  const int ncol = moments ? 8 : 4;
  Rcpp::NumericMatrix out(np, ncol);
  double* o = out.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (long k = 0; k < (long)tasks.size(); ++k) {
    const Task& tk = tasks[k];
    MatD Uab = A.middleCols(i0 + tk.ib, tk.nb).transpose() * A.middleCols(tk.j0, tk.nj);
    MatD T1 = CT.middleRows(i0 + tk.ib, tk.nb) * Cd.middleRows(tk.j0, tk.nj).transpose();
    MatD T2 = S[2].middleCols(tk.ib, tk.nb).transpose() * K[2].middleCols(tk.j0, tk.nj);
    MatD T3 = S[3].middleCols(tk.ib, tk.nb).transpose() * K[3].middleCols(tk.j0, tk.nj);
    MatD T4 = S[4].middleCols(tk.ib, tk.nb).transpose() * K[4].middleCols(tk.j0, tk.nj);
    for (int jj = 0; jj < tk.nj; ++jj) {
      const int jg = tk.j0 + jj;
      for (int ii = 0; ii < tk.nb; ++ii) {
        const int il = tk.ib + ii, ig = i0 + il;
        if (jg <= ig) continue;
        write_pair(o, np, offset[il] + (jg - ig - 1), moments, Uab(ii, jj),
                   T1(ii, jj), T2(ii, jj), T3(ii, jj), T4(ii, jj));
      }
    }
  }
  Rcpp::colnames(out) = pair_columns(moments);
  return out;
}
