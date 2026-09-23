#define EIGEN_DONT_PARALLELIZE
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <climits>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::export]]
arma::mat mgcvst_pair_trace_powers_cpp(const Rcpp::List& matrixList,
                                       const Rcpp::IntegerMatrix& pairs,
                                       int maxPower = 4,
                                       int threads = 1) {
  int n = matrixList.size();
  if (n < 1) {
    Rcpp::stop("matrixList must contain at least one matrix.");
  }
  if (pairs.ncol() != 2 || pairs.nrow() < 1) {
    Rcpp::stop("pairs must be a non-empty two-column integer matrix.");
  }
  if (maxPower < 1 || maxPower > 4) {
    Rcpp::stop("maxPower must be an integer from 1 to 4.");
  }
  if (threads < 1) {
    Rcpp::stop("threads must be a positive integer.");
  }
#ifndef _OPENMP
  if (threads > 1) {
    Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
  }
#endif

  Rcpp::NumericMatrix first(matrixList[0]);
  int q = first.nrow();
  if (q < 1 || first.ncol() != q) {
    Rcpp::stop("Every matrix must be non-empty and square.");
  }
  std::vector<double*> matrixPointers(n);
  matrixPointers[0] = first.begin();
  for (int i = 1; i < n; ++i) {
    Rcpp::NumericMatrix current(matrixList[i]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every matrix must have the same square dimension.");
    }
    matrixPointers[i] = current.begin();
  }
  for (int k = 0; k < pairs.nrow(); ++k) {
    int i = pairs(k, 0) - 1;
    int j = pairs(k, 1) - 1;
    if (i < 0 || i >= n || j < 0 || j >= n) {
      Rcpp::stop("pairs contains an index outside matrixList.");
    }
  }

  arma::mat out(pairs.nrow(), maxPower, arma::fill::none);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int k = 0; k < pairs.nrow(); ++k) {
    int i = pairs(k, 0) - 1;
    int j = pairs(k, 1) - 1;
    // Keep BLAS calls out of the outer OpenMP region. The conda pthread
    // OpenBLAS build cannot safely create its own workers here.
    Eigen::Map<const Eigen::MatrixXd> left(matrixPointers[i], q, q);
    Eigen::Map<const Eigen::MatrixXd> right(matrixPointers[j], q, q);
    Eigen::MatrixXd product = left * right;
    out(k, 0) = product.trace();

    if (maxPower >= 2) {
      out(k, 1) = (product.array() * product.transpose().array()).sum();
      if (maxPower >= 3) {
        Eigen::MatrixXd product2 = product * product;
        out(k, 2) =
          (product2.array() * product.transpose().array()).sum();
        if (maxPower >= 4) {
          out(k, 3) =
            (product2.array() * product2.transpose().array()).sum();
        }
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_pair_lowrank_cpp(
    const Rcpp::NumericMatrix& A, const Rcpp::List& left,
    const Rcpp::List& right, const Rcpp::NumericMatrix& scale,
    const Rcpp::IntegerMatrix& pairs, int threads = 1) {
  const int genes = A.ncol(), q = A.nrow();
  if (genes < 1 || q < 1 || left.size() != 4 || right.size() != 4 ||
      scale.nrow() != genes || scale.ncol() != 4 || pairs.ncol() != 2 ||
      threads < 1) Rcpp::stop("Incompatible low-rank pair inputs.");
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
#endif
  std::vector<const double*> lp(4), rp(4);
  std::vector<int> ranks(4);
  for (int k = 0; k < 4; ++k) {
    Rcpp::NumericMatrix L(left[k]), R(right[k]);
    if (L.ncol() != genes || R.ncol() != genes ||
        L.nrow() != R.nrow() || L.nrow() < 1)
      Rcpp::stop("Low-rank factors must align with the score columns.");
    lp[k] = L.begin(); rp[k] = R.begin(); ranks[k] = L.nrow();
  }
  for (int k = 0; k < pairs.nrow(); ++k) {
    if (pairs(k, 0) < 1 || pairs(k, 0) > genes ||
        pairs(k, 1) < 1 || pairs(k, 1) > genes)
      Rcpp::stop("pairs contains an index outside the score columns.");
  }
  const int n = pairs.nrow();
  const double* a = A.begin();
  const double* s = scale.begin();
  const int* pair = pairs.begin();
  Rcpp::NumericMatrix out(n, 5);
  double* output = out.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int row = 0; row < n; ++row) {
    const int i = pair[row] - 1, j = pair[row + n] - 1;
    double score = 0;
    for (int d = 0; d < q; ++d)
      score += a[static_cast<size_t>(i) * q + d] * a[static_cast<size_t>(j) * q + d];
    output[row] = score;
    for (int k = 0; k < 4; ++k) {
      double moment = 0;
      const double* l = lp[k] + static_cast<size_t>(i) * ranks[k];
      const double* r = rp[k] + static_cast<size_t>(j) * ranks[k];
      for (int d = 0; d < ranks[k]; ++d) moment += l[d] * r[d];
      output[row + static_cast<size_t>(k + 1) * n] =
        moment * s[i + static_cast<size_t>(k) * genes] * s[j + static_cast<size_t>(k) * genes];
    }
  }
  return out;
}

// Rows are gene-major: row g * referenceList.size() + r contains
// the four trace powers for matrixList[g] and referenceList[r].
// [[Rcpp::export]]
Rcpp::NumericMatrix mgcvst_landmark_trace_cpp(
    const Rcpp::List& matrixList, const Rcpp::List& referenceList,
    int threads = 1, bool float32 = true) {
  const int genes = matrixList.size();
  const int references = referenceList.size();
  if (genes < 1 || references < 1) {
    Rcpp::stop("matrixList and referenceList must be non-empty.");
  }
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
#ifndef _OPENMP
  if (threads > 1) {
    Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
  }
#endif
  if (genes > INT_MAX / references) {
    Rcpp::stop("The landmark trace result has too many rows.");
  }

  Rcpp::NumericMatrix first(matrixList[0]);
  const int q = first.nrow();
  if (q < 1 || first.ncol() != q) {
    Rcpp::stop("Every input matrix must be non-empty and square.");
  }
  std::vector<const double*> genePointers(genes), referencePointers(references);
  genePointers[0] = first.begin();
  for (int g = 1; g < genes; ++g) {
    Rcpp::NumericMatrix current(matrixList[g]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every input matrix must have the same square dimension.");
    }
    genePointers[g] = current.begin();
  }
  for (int r = 0; r < references; ++r) {
    Rcpp::NumericMatrix current(referenceList[r]);
    if (current.nrow() != q || current.ncol() != q) {
      Rcpp::stop("Every input matrix must have the same square dimension.");
    }
    referencePointers[r] = current.begin();
  }

  const int rows = genes * references;
  Rcpp::NumericMatrix out(rows, 4);
  double* output = out.begin();

  if (float32) {
    std::vector<Eigen::MatrixXf> refFloat;
    refFloat.reserve(references);
    for (int r = 0; r < references; ++r) {
      Eigen::Map<const Eigen::MatrixXd> ref(referencePointers[r], q, q);
      refFloat.emplace_back(ref.cast<float>());
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(threads < genes ? threads : genes)
#endif
    {
      Eigen::MatrixXf gene(q, q), product(q, q), product2(q, q);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int g = 0; g < genes; ++g) {
        Eigen::Map<const Eigen::MatrixXd> source(genePointers[g], q, q);
        gene = source.cast<float>();
        for (int r = 0; r < references; ++r) {
          product.noalias() = gene * refFloat[r];
          product2.noalias() = product * product;
          double t1 = 0, t2 = 0, t3 = 0, t4 = 0;
          for (int col = 0; col < q; ++col) {
            t1 += static_cast<double>(product(col, col));
            for (int row = 0; row < q; ++row) {
              const double p = product(row, col);
              const double pt = product(col, row);
              const double p2 = product2(row, col);
              t2 += p * pt;
              t3 += p2 * pt;
              t4 += p2 * static_cast<double>(product2(col, row));
            }
          }
          const int index = g * references + r;
          output[index] = t1;
          output[index + rows] = t2;
          output[index + 2 * static_cast<size_t>(rows)] = t3;
          output[index + 3 * static_cast<size_t>(rows)] = t4;
        }
      }
    }
  } else {
#ifdef _OPENMP
#pragma omp parallel num_threads(threads < genes ? threads : genes)
#endif
    {
      Eigen::MatrixXd product(q, q), product2(q, q);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int g = 0; g < genes; ++g) {
        Eigen::Map<const Eigen::MatrixXd> gene(genePointers[g], q, q);
        for (int r = 0; r < references; ++r) {
          Eigen::Map<const Eigen::MatrixXd> ref(referencePointers[r], q, q);
          product.noalias() = gene * ref;
          product2.noalias() = product * product;
          double t1 = 0, t2 = 0, t3 = 0, t4 = 0;
          for (int col = 0; col < q; ++col) {
            t1 += product(col, col);
            for (int row = 0; row < q; ++row) {
              const double p = product(row, col);
              const double pt = product(col, row);
              const double p2 = product2(row, col);
              t2 += p * pt;
              t3 += p2 * pt;
              t4 += p2 * product2(col, row);
            }
          }
          const int index = g * references + r;
          output[index] = t1;
          output[index + rows] = t2;
          output[index + 2 * static_cast<size_t>(rows)] = t3;
          output[index + 3 * static_cast<size_t>(rows)] = t4;
        }
      }
    }
  }
  return out;
}
