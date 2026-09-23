# Read-only system probes. Missing signals remain unknown, not zero.
.mgcvst_memory_lines <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path) ||
      file.access(path, 4L) != 0L) return(character())
  readLines(path, warn = FALSE)
}

.mgcvst_memory_number <- function(x) {
  if (length(x) != 1L || is.na(x) || !grepl("^[0-9]+([.][0-9]+)?$", x)) {
    return(NA_real_)
  }
  z <- as.numeric(x)
  if (is.finite(z)) z else NA_real_
}

.mgcvst_memory_kb <- function(lines, key) {
  z <- lines[startsWith(lines, paste0(key, ":"))]
  if (length(z) != 1L) return(NA_real_)
  fields <- strsplit(trimws(sub("^[^:]+:", "", z)), "[[:space:]]+")[[1L]]
  if (length(fields) != 2L || fields[2L] != "kB") return(NA_real_)
  .mgcvst_memory_number(fields[1L]) * 1024
}

.mgcvst_memory_unescape_mount <- function(x) {
  for (code in c("040", "011", "012", "134")) {
    x <- gsub(paste0("\\", code), intToUtf8(strtoi(code, base = 8L)),
              x, fixed = TRUE)
  }
  x
}

.mgcvst_slurm_task_count <- function(env) {
  for (key in c("SLURM_NTASKS_PER_NODE", "SLURM_TASKS_PER_NODE")) {
    x <- unname(env[key])
    if (!length(x) || is.na(x) || !nzchar(x)) next
    token <- strsplit(x, ",", fixed = TRUE)[[1L]]
    good <- grepl("^[0-9]+(\\(x[0-9]+\\))?$", token)
    if (all(good)) {
      n <- as.numeric(sub("\\(.*$", "", token))
      if (all(is.finite(n)) && all(n > 0)) return(max(n))
    }
  }
  ntasks <- .mgcvst_memory_number(unname(env["SLURM_NTASKS"]))
  nodes <- .mgcvst_memory_number(unname(env["SLURM_JOB_NUM_NODES"]))
  if (is.finite(ntasks) && ntasks > 0 && is.finite(nodes) && nodes > 0) {
    return(ceiling(ntasks / nodes))
  }
  NA_real_
}

# Resolve the process cgroup against its mount, including constrained ancestors.
.mgcvst_cgroup_headroom <- function(cgroups, mounts) {
  out <- numeric()
  for (line in mounts) {
    fields <- strsplit(line, " - ", fixed = TRUE)[[1L]]
    if (length(fields) != 2L) next
    left <- strsplit(fields[1L], " ", fixed = TRUE)[[1L]]
    right <- strsplit(fields[2L], " ", fixed = TRUE)[[1L]]
    if (length(left) < 5L || length(right) < 3L) next
    v2 <- identical(right[1L], "cgroup2")
    v1 <- identical(right[1L], "cgroup") &&
      "memory" %in% strsplit(right[3L], ",", fixed = TRUE)[[1L]]
    if (!v1 && !v2) next
    root <- .mgcvst_memory_unescape_mount(left[4L])
    mount <- sub("/+$", "", .mgcvst_memory_unescape_mount(left[5L]))
    if (!nzchar(mount)) mount <- "/"
    for (group in cgroups) {
      z <- strsplit(group, ":", fixed = TRUE)[[1L]]
      if (length(z) < 3L ||
          (v2 && nzchar(z[2L])) ||
          (v1 && !("memory" %in% strsplit(z[2L], ",", fixed = TRUE)[[1L]]))) next
      path <- paste(z[-c(1L, 2L)], collapse = ":")
      if (identical(path, root)) {
        relative <- ""
      } else if (identical(root, "/")) {
        relative <- path
      } else if (startsWith(path, paste0(root, "/"))) {
        relative <- substring(path, nchar(root) + 1L)
      } else if (identical(path, "/")) {
        relative <- "" # cgroup namespace exposes its mounted root as '/'.
      } else next
      if (any(strsplit(relative, "/", fixed = TRUE)[[1L]] == "..")) next
      current <- paste0(mount, relative)
      repeat {
        limit <- .mgcvst_memory_number(.mgcvst_memory_lines(file.path(
          current, if (v2) "memory.max" else "memory.limit_in_bytes")))
        usage <- .mgcvst_memory_number(.mgcvst_memory_lines(file.path(
          current, if (v2) "memory.current" else "memory.usage_in_bytes")))
        if (is.finite(limit) && limit > 0 && limit < 2^60 &&
            is.finite(usage) && usage >= 0) {
          out[paste0(if (v2) "cgroup2:" else "cgroup1:", current)] <-
            max(0, limit - usage)
        }
        if (identical(current, mount)) break
        parent <- dirname(current)
        if (identical(parent, current) ||
            !(identical(parent, mount) || startsWith(parent, paste0(mount, "/")))) break
        current <- parent
      }
    }
  }
  out
}

.mgcvst_memory_probe <- function(proc = "/proc", env = Sys.getenv(),
                                  native = mgcvst_memory_status_cpp()) {
  signals <- numeric()
  resident <- unname(native["resident"])
  if (is.finite(native["available"])) signals["system_available"] <- native["available"]
  available <- .mgcvst_memory_kb(.mgcvst_memory_lines(file.path(proc, "meminfo")), "MemAvailable")
  rss <- .mgcvst_memory_kb(.mgcvst_memory_lines(file.path(proc, "self/status")), "VmRSS")
  if (is.finite(available)) signals["MemAvailable"] <- available
  if (is.finite(rss)) resident <- rss
  if (!is.finite(resident) || resident < 0) resident <- NA_real_
  signals <- c(signals, .mgcvst_cgroup_headroom(
    .mgcvst_memory_lines(file.path(proc, "self/cgroup")),
    .mgcvst_memory_lines(file.path(proc, "self/mountinfo"))))

  limit <- .mgcvst_memory_number(unname(env["SLURM_MEM_PER_NODE"])) * 1024^2
  if (!is.finite(limit)) {
    per_cpu <- .mgcvst_memory_number(unname(env["SLURM_MEM_PER_CPU"])) * 1024^2
    cpus <- .mgcvst_memory_number(unname(env["SLURM_CPUS_ON_NODE"]))
    if (!is.finite(cpus)) {
      cpu_list <- unname(env["SLURM_JOB_CPUS_PER_NODE"])
      if (length(cpu_list) == 1L && !is.na(cpu_list) && nzchar(cpu_list)) {
        token <- strsplit(cpu_list, ",", fixed = TRUE)[[1L]]
        good <- grepl("^[0-9]+(\\(x[0-9]+\\))?$", token)
        if (all(good)) cpus <- max(as.numeric(sub("\\(.*$", "", token)))
      }
    }
    if (is.finite(cpus) && cpus > 0) limit <- per_cpu * cpus
  }
  # Share node memory conservatively when the job has multiple local tasks.
  tasks <- .mgcvst_slurm_task_count(env)
  if (is.finite(limit) && limit > 0 && is.finite(resident)) {
    if (is.finite(tasks) && tasks > 0) limit <- limit / tasks
    signals["slurm_allocation_remaining"] <- max(0, limit - resident)
  }
  signals <- signals[is.finite(signals) & signals >= 0]
  list(available = if (length(signals)) min(signals) else NA_real_,
       resident = resident, signals = signals,
       source = if (length(signals)) names(signals)[which.min(signals)] else "unknown")
}

# Budget is additional pair-stage memory, after fit/basis/result allocation.
.mgcvst_inla_memory_plan <- function(fit, basis, pairs, threads,
                                      memory_gb = NULL, resident_cache = 0,
                                      state_bytes = NULL, probe = .mgcvst_memory_probe()) {
  if (!is.null(memory_gb) && (!is.numeric(memory_gb) || length(memory_gb) != 1L ||
      !is.finite(memory_gb) || memory_gb <= 0)) {
    stop("memory_gb must be NULL or one positive finite number of GiB.")
  }
  q <- ncol(fit$score_sparse$Q)
  r <- basis$rank
  p <- ncol(fit$geometry$nuisance_design)
  if (is.null(p)) p <- 0L
  n <- nrow(fit$working_variance)
  if (is.null(n)) n <- 0L
  estimate <- 8 * (r^2 + r + p^2) + 2048
  if (!is.null(state_bytes)) estimate <- max(estimate, state_bytes)
  pair_work <- 3 * 8 * r^2 * min(threads, pairs)
  prepare_work <- 8 * (6 * q^2 + 4 * q * p + 4 * p^2 + 2 * n + 4 * r^2)
  # Account for retained pair result chunks and final combination, not one R frame per pair.
  reserve <- max(pair_work, prepare_work) + 256 * pairs + 65536
  detected <- if (is.finite(probe$available)) probe$available + resident_cache else NA_real_
  budget <- if (is.null(memory_gb)) {
    if (is.finite(detected)) 0.8 * detected else 512 * 1024^2
  } else {
    if (is.finite(detected)) min(memory_gb * 1024^3, 0.8 * detected) else memory_gb * 1024^3
  }
  capacity <- floor((budget - reserve) / estimate)
  if (!is.finite(capacity) || capacity < 2) {
    stop("Insufficient detected or requested memory for two INLA score states and the working buffers. ",
         "Reduce threads or increase the job memory allocation.")
  }
  list(budget_bytes = budget, cache_bytes = budget - reserve,
       reserve_bytes = reserve, state_bytes = estimate,
       capacity = min(capacity, .Machine$integer.max), requested_gb = memory_gb,
       available_bytes = probe$available, source = probe$source,
       signals = probe$signals, fallback = !is.finite(detected))
}
