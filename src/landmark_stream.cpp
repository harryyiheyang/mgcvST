#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include "score_state_io.h"
#include <climits>
#include <numeric>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
void mgcvst_state_write_cpp(const std::string& path,
                            const std::string& signature,
                            const std::string& feature_id,
                            const Rcpp::NumericVector& a,
                            const Rcpp::NumericMatrix& M,
                            SEXP width,
                            const std::string& error = "") {
  mgcvst_io::State state;
  state.error = error;
  if (error.empty()) {
    const int q = M.nrow();
    if (q < 1 || M.ncol() != q || a.size() != q) {
      Rcpp::stop("Native state needs a square M and aligned score vector.");
    }
    state.q = static_cast<uint64_t>(q);
    state.a.assign(a.begin(), a.end());
    if (TYPEOF(width) == REALSXP) {
      Rcpp::NumericVector values(width);
      state.width.assign(values.begin(), values.end());
    } else if (TYPEOF(width) == INTSXP) {
      Rcpp::IntegerVector values(width);
      for (int k = 0; k < values.size(); ++k) {
        if (values[k] == NA_INTEGER) Rcpp::stop("Native widths cannot be NA.");
        state.width.push_back(values[k]);
      }
    } else {
      Rcpp::stop("Native widths must be numeric or integer.");
    }
    SEXP name_attr = Rf_getAttrib(width, R_NamesSymbol);
    if (name_attr == R_NilValue) {
      state.width_names.assign(state.width.size(), "");
    } else {
      Rcpp::CharacterVector names(name_attr);
      if (names.size() != static_cast<int>(state.width.size())) Rcpp::stop("Invalid native width names.");
      for (int k = 0; k < names.size(); ++k) {
        if (names[k] == NA_STRING) Rcpp::stop("Native width names cannot be NA.");
        state.width_names.push_back(Rcpp::as<std::string>(names[k]));
      }
    }
    if (!std::all_of(state.a.begin(), state.a.end(), [](double z) { return std::isfinite(z); }) ||
        !std::all_of(state.width.begin(), state.width.end(), [](double z) { return std::isfinite(z); })) {
      Rcpp::stop("Native state has non-finite scores or widths.");
    }
    state.upper.reserve(static_cast<size_t>(q) * (static_cast<size_t>(q) + 1) / 2);
    for (int col = 0; col < q; ++col) {
      for (int row = 0; row <= col; ++row) {
        const double upper = M(row, col), lower = M(col, row);
        if (!std::isfinite(upper) || !std::isfinite(lower) ||
            std::abs(upper - lower) > 1e-8 * std::max({1.0, std::abs(upper), std::abs(lower)})) {
          Rcpp::stop("Native score matrix must be finite and symmetric.");
        }
        state.upper.push_back(upper);
      }
    }
  }
  mgcvst_io::write_state(path, signature, feature_id, state);
}

// [[Rcpp::export]]
Rcpp::List mgcvst_state_read_cpp(const std::string& path,
                                const std::string& signature,
                                const std::string& feature_id) {
  mgcvst_io::State state = mgcvst_io::read_state(path, signature, feature_id);
  if (!state.error.empty()) return Rcpp::List::create(Rcpp::_["error"] = state.error);
  if (state.q > INT_MAX) Rcpp::stop("Native score matrix exceeds R dimensions.");
  const int q = static_cast<int>(state.q);
  Rcpp::NumericMatrix M(q, q);
  size_t at = 0;
  for (int col = 0; col < q; ++col) {
    for (int row = 0; row <= col; ++row) {
      const double value = state.upper[at++];
      M(row, col) = value; M(col, row) = value;
    }
  }
  Rcpp::NumericVector width(state.width.begin(), state.width.end());
  bool named = std::any_of(state.width_names.begin(), state.width_names.end(),
    [](const std::string& name) { return !name.empty(); });
  if (named) {
    Rcpp::CharacterVector names(state.width_names.size());
    for (size_t k = 0; k < state.width_names.size(); ++k) names[k] = state.width_names[k];
    width.attr("names") = names;
  }
  return Rcpp::List::create(Rcpp::_["a"] = state.a, Rcpp::_["M"] = M,
                            Rcpp::_["width"] = width);
}

// Gene-major input queues are processed by dynamic OpenMP scheduling. Each
// worker holds one packed double shard and three float matrices at most.
// [[Rcpp::export]]
Rcpp::List mgcvst_landmark_stream_cpp(
    const Rcpp::CharacterVector& paths,
    const Rcpp::CharacterVector& feature_ids,
    const std::string& state_signature,
    const Rcpp::List& references,
    const Rcpp::CharacterVector& output_paths,
    const std::string& summary_signature,
    int threads = 1) {
  const int genes = paths.size(), n_ref = references.size();
  if (genes < 1 || n_ref < 1 || feature_ids.size() != genes ||
      output_paths.size() != genes || threads < 1 ||
      state_signature.empty() || summary_signature.empty()) {
    Rcpp::stop("Invalid native landmark stream queue or signatures.");
  }
#ifndef _OPENMP
  if (threads > 1) Rcpp::stop("mgcvST was compiled without OpenMP support; use threads = 1.");
#endif
  std::vector<std::string> input(genes), ids(genes), output(genes);
  std::unordered_set<std::string> output_seen;
  for (int g = 0; g < genes; ++g) {
    if (STRING_ELT(paths, g) == NA_STRING || STRING_ELT(feature_ids, g) == NA_STRING ||
        STRING_ELT(output_paths, g) == NA_STRING) {
      Rcpp::stop("Native landmark paths and IDs cannot be NA.");
    }
    input[g] = Rcpp::as<std::string>(paths[g]);
    ids[g] = Rcpp::as<std::string>(feature_ids[g]);
    output[g] = Rcpp::as<std::string>(output_paths[g]);
    if (input[g].empty() || ids[g].empty() || output[g].empty() ||
        input[g] == output[g]) Rcpp::stop("Invalid native landmark path or feature ID.");
    if (!output_seen.insert(output[g]).second) {
      Rcpp::stop("Native landmark output paths must be distinct.");
    }
  }
  Rcpp::NumericMatrix first(references[0]);
  const int q = first.nrow();
  if (q < 1 || first.ncol() != q) Rcpp::stop("References must be non-empty square matrices.");
  std::vector<Eigen::MatrixXf> ref_float;
  ref_float.reserve(n_ref);
  for (int r = 0; r < n_ref; ++r) {
    Rcpp::NumericMatrix matrix(references[r]);
    if (matrix.nrow() != q || matrix.ncol() != q) {
      Rcpp::stop("References must have one square dimension.");
    }
    Eigen::Map<const Eigen::MatrixXd> source(matrix.begin(), q, q);
    if (!source.allFinite()) Rcpp::stop("Reference matrices must be finite.");
    ref_float.emplace_back(source.cast<float>());
  }
  std::vector<int> completed(genes, 0), source_errors(genes, 0);
  std::vector<std::string> errors(genes), source_messages(genes);
  const int workers = std::min(threads, genes);
#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
  {
    Eigen::MatrixXf gene, product, product2;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int g = 0; g < genes; ++g) {
      try {
        // Keep allocation failures inside the worker's exception boundary.
        gene.resize(q, q); product.resize(q, q); product2.resize(q, q);
        mgcvst_io::State state = mgcvst_io::read_state(input[g], state_signature, ids[g]);
        mgcvst_io::Trace trace;
        trace.n_ref = static_cast<uint64_t>(n_ref);
        if (!state.error.empty()) {
          trace.error = state.error;
          source_errors[g] = 1;
          source_messages[g] = state.error;
        } else {
          if (state.q != static_cast<uint64_t>(q) || state.a.size() != static_cast<size_t>(q)) {
            throw std::runtime_error("Native score state and references have different dimensions.");
          }
          trace.a = std::move(state.a);
          size_t at = 0;
          for (int col = 0; col < q; ++col) {
            for (int row = 0; row <= col; ++row) {
              const float value = static_cast<float>(state.upper[at++]);
              if (!std::isfinite(value)) throw std::runtime_error("Score matrix overflows float32.");
              gene(row, col) = value; gene(col, row) = value;
            }
          }
          std::vector<double>().swap(state.upper);
          trace.cross.resize(static_cast<size_t>(n_ref) * 4);
          for (int r = 0; r < n_ref; ++r) {
            product.noalias() = gene * ref_float[r];
            product2.noalias() = product * product;
            double t1 = 0, t2 = 0, t3 = 0, t4 = 0;
            for (int col = 0; col < q; ++col) {
              t1 += static_cast<double>(product(col, col));
              for (int row = 0; row < q; ++row) {
                const double p = product(row, col), pt = product(col, row);
                const double p2 = product2(row, col);
                t2 += p * pt;
                t3 += p2 * pt;
                t4 += p2 * static_cast<double>(product2(col, row));
              }
            }
            const size_t k = static_cast<size_t>(r);
            trace.cross[k] = t1;
            trace.cross[k + n_ref] = t2;
            trace.cross[k + 2 * static_cast<size_t>(n_ref)] = t3;
            trace.cross[k + 3 * static_cast<size_t>(n_ref)] = t4;
          }
          if (!std::all_of(trace.cross.begin(), trace.cross.end(),
                [](double z) { return std::isfinite(z); })) {
            throw std::runtime_error("Native landmark trace is non-finite.");
          }
        }
        mgcvst_io::write_trace(output[g], summary_signature, ids[g], trace);
        completed[g] = 1;
      } catch (const std::exception& e) {
        errors[g] = e.what();
      }
    }
  }
  const int done = std::count(completed.begin(), completed.end(), 1);
  if (done != genes) {
    for (int g = 0; g < genes; ++g) {
      if (!completed[g]) {
        Rcpp::stop("Native landmark stream completed %d of %d genes; %s: %s",
                   done, genes, ids[g].c_str(), errors[g].c_str());
      }
    }
  }
  const double qq = static_cast<double>(q) * q;
  // Includes packed-double read buffers, three float matrices, and allowance
  // for Eigen GEMM packing; R's original double references are additional.
  const double workspace = 4 * qq * n_ref + 24 * qq * workers +
    8 * static_cast<double>(q) * workers + 32 * static_cast<double>(n_ref) * workers;
  Rcpp::CharacterVector source_error_message(genes);
  for (int g = 0; g < genes; ++g) {
    source_error_message[g] = source_errors[g] ?
      Rcpp::String(source_messages[g]) : Rcpp::String(NA_STRING);
  }
  return Rcpp::List::create(Rcpp::_["completed"] = done,
                            Rcpp::_["source_errors"] = std::accumulate(source_errors.begin(), source_errors.end(), 0),
                            Rcpp::_["source_error_message"] = source_error_message,
                            Rcpp::_["workspace_bytes_estimate"] = workspace,
                            Rcpp::_["threads"] = workers);
}

// [[Rcpp::export]]
Rcpp::List mgcvst_trace_read_cpp(const std::string& path,
                                const std::string& signature,
                                const std::string& feature_id,
                                int n_ref) {
  if (n_ref < 1) Rcpp::stop("n_ref must be positive.");
  mgcvst_io::Trace trace = mgcvst_io::read_trace(path, signature, feature_id,
                                                  static_cast<uint64_t>(n_ref));
  if (!trace.error.empty()) return Rcpp::List::create(Rcpp::_["error"] = trace.error);
  Rcpp::NumericMatrix cross(n_ref, 4);
  std::copy(trace.cross.begin(), trace.cross.end(), cross.begin());
  return Rcpp::List::create(Rcpp::_["a"] = trace.a,
                            Rcpp::_["cross"] = cross,
                            Rcpp::_["self"] = R_NilValue);
}
