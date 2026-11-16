// Run-time choice of the SIMD kernels in la_kernel.cpp (the dense products and
// symmetric eigenvalue problems of the pair tests; compiled for generic,
// AVX2+FMA and AVX-512 by src/Makevars). la() returns the best table the CPU
// supports, chosen once per process; the environment variable
// MGCVST_SIMD=avx512|avx2|generic caps the choice. The generic copy is compiled
// with the package's ordinary flags. This is the package's only SIMD dispatch.
#ifndef MGCVST_LA_DISPATCH_H
#define MGCVST_LA_DISPATCH_H

#include <cstdlib>
#include <string>

namespace mgcvst {

struct LaFns {
  // Work space of one thread: q x q matrices for the trace moments when
  // `moments` is set, k x k matrices and the eigenvalue solver.
  void* (*work_new)(int q, int k, int moments);
  void (*work_free)(void* work);
  // Left matrix of the trace moments: Ln = sym(H) / scale (q x q), kept.
  void (*set_left)(void* work, const double* H, double scale);
  // v[s - 1] = tr((Ln Rn)^s), s = 1..4, Rn = sym(H) / scale.
  void (*moments)(void* work, const double* H, double scale, double* v);
  // Singular values (descending, k) of Gl' Gr, Gl and Gr q x k: square roots of
  // the eigenvalues of the Gram matrix. 0 done, 1 non-finite product, 2 the
  // eigenvalue solver failed.
  int (*spectrum_tn)(void* work, const double* Gl, const double* Gr, int q, double* sv);
  // The same for Ri Rj', Ri and Rj k x k.
  int (*spectrum_nt)(void* work, const double* Ri, const double* Rj, double* sv);
  // Out (M x N) = A (M x K) B (K x N), column-major doubles, one serial call.
  void (*dgemm_nn)(const double* A, int M, int K, const double* B, int N, double* Out);
  // Out (M x N) = A' B, A K x M and B K x N, column-major floats, one serial call.
  void (*sgemm_tn)(const float* A, long K, int M, const float* B, int N, float* Out);
  const char* isa;
};

}  // namespace mgcvst

extern const mgcvst::LaFns mgcvst_la_generic;
#if defined(MGCVST_LA_X86)
extern const mgcvst::LaFns mgcvst_la_avx2;
extern const mgcvst::LaFns mgcvst_la_avx512;
#endif

namespace mgcvst {

inline const LaFns& pick_la() {
  const char* env = std::getenv("MGCVST_SIMD");
  const std::string want = env ? env : "";
#if defined(MGCVST_LA_X86) && (defined(__GNUC__) || defined(__clang__))
  __builtin_cpu_init();
  const bool a512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
                    __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512bw") &&
                    __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
  const bool a2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
  if (a512 && (want.empty() || want == "avx512")) return mgcvst_la_avx512;
  if (a2 && (want.empty() || want == "avx512" || want == "avx2")) return mgcvst_la_avx2;
#endif
  return mgcvst_la_generic;
}
inline const LaFns& la() { static const LaFns& t = pick_la(); return t; }

// Work space of one thread, released with the thread's block.
struct LaWork {
  const LaFns& f;
  void* p;
  LaWork(const LaFns& fns, int q, int k, bool moments)
    : f(fns), p(fns.work_new(q, k, moments ? 1 : 0)) {}
  ~LaWork() { f.work_free(p); }
  LaWork(const LaWork&) = delete;
  LaWork& operator=(const LaWork&) = delete;
};

}  // namespace mgcvst
#endif
