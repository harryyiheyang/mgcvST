memory_fixture_write <- function(root, path, text) {
  target <- file.path(root, path)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, target, useBytes = TRUE)
  target
}

memory_proc_fixture <- function(root, available_kb = NULL, rss_kb = NULL,
                                cgroup = character(), mountinfo = character()) {
  if (!is.null(available_kb)) {
    memory_fixture_write(root, "meminfo",
                         paste0("MemAvailable: ",
                                format(available_kb, scientific = FALSE, trim = TRUE),
                                " kB"))
  }
  if (!is.null(rss_kb)) {
    memory_fixture_write(root, "self/status",
                         paste0("VmRSS: ", format(rss_kb, scientific = FALSE,
                                                   trim = TRUE), " kB"))
  }
  if (length(cgroup)) memory_fixture_write(root, "self/cgroup", cgroup)
  if (length(mountinfo)) memory_fixture_write(root, "self/mountinfo", mountinfo)
  root
}

test_that("memory parsers preserve unknown values", {
  expect_true(is.na(mgcvST:::.mgcvst_memory_number("max")))
  expect_true(is.na(mgcvST:::.mgcvst_memory_number("-1")))
  expect_true(is.na(mgcvST:::.mgcvst_memory_kb("MemTotal: 4 kB", "MemAvailable")))
  expect_equal(mgcvST:::.mgcvst_memory_unescape_mount("/a\\040b\\011c"),
               paste0("/a b", "\tc"))
})

test_that("cgroup v2 uses the tightest nested remaining limit", {
  root <- tempfile("memory-v2-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  memory_fixture_write(root, "cg/job/step/memory.max", "100000000000")
  memory_fixture_write(root, "cg/job/step/memory.current", "80000000000")
  memory_fixture_write(root, "cg/job/memory.max", "15000000000")
  memory_fixture_write(root, "cg/job/memory.current", "12000000000")
  memory_fixture_write(root, "cg/memory.max", "max")
  memory_fixture_write(root, "cg/memory.current", "100")
  got <- mgcvST:::.mgcvst_cgroup_headroom(
    "0::/job/step",
    paste0("29 23 0:26 / ", normalizePath(file.path(root, "cg"), winslash = "/"),
           " rw - cgroup2 cgroup rw")
  )
  expect_identical(unname(min(got)), 3000000000)
  memory_fixture_write(root, "cg/job/memory.current", "15000000000")
  got <- mgcvST:::.mgcvst_cgroup_headroom(
    "0::/job/step",
    paste0("29 23 0:26 / ", normalizePath(file.path(root, "cg"), winslash = "/"),
           " rw - cgroup2 cgroup rw")
  )
  expect_identical(unname(min(got)), 0)
})

test_that("cgroup v1 handles unlimited sentinel and namespace mounts", {
  root <- tempfile("memory v1 ")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  mount <- file.path(root, "memory")
  memory_fixture_write(root, "memory/memory.limit_in_bytes", "9223372036854771712")
  memory_fixture_write(root, "memory/memory.usage_in_bytes", "100")
  got <- mgcvST:::.mgcvst_cgroup_headroom(
    "5:cpu,memory:/",
    paste0("31 23 0:27 /tenant/job ", gsub(" ", "\\040", normalizePath(mount, winslash = "/"), fixed = TRUE),
           " rw - cgroup cgroup rw,memory")
  )
  expect_length(got, 0L)
  memory_fixture_write(root, "memory/memory.limit_in_bytes", "2000000000")
  memory_fixture_write(root, "memory/memory.usage_in_bytes", "1250000000")
  got <- mgcvST:::.mgcvst_cgroup_headroom(
    "5:cpu,memory:/tenant/job",
    paste0("31 23 0:27 /tenant/job ", gsub(" ", "\\040", normalizePath(mount, winslash = "/"), fixed = TRUE),
           " rw - cgroup cgroup rw,memory")
  )
  expect_identical(unname(got), 750000000)
})

test_that("Slurm memory uses node CPU and task allocation conservatively", {
  root <- tempfile("memory-slurm-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  proc <- memory_proc_fixture(root, available_kb = 200000000,
                              rss_kb = 4000000)
  expect_true(file.exists(file.path(proc, "self/status")))
  expect_identical(readLines(file.path(proc, "self/status")),
                   "VmRSS: 4000000 kB")
  expect_equal(mgcvST:::.mgcvst_memory_kb(
    mgcvST:::.mgcvst_memory_lines(file.path(proc, "self/status")), "VmRSS"
  ), 4096000000)
  native <- c(available = NA_real_, resident = NA_real_)
  per_node <- mgcvST:::.mgcvst_memory_probe(
    proc, c(SLURM_MEM_PER_NODE = "150000", SLURM_NTASKS_PER_NODE = "2(x3)"), native
  )
  expect_equal(unname(per_node$signals["slurm_allocation_remaining"]),
               74547200000)
  expect_equal(per_node$available, 74547200000)
  per_cpu <- mgcvST:::.mgcvst_memory_probe(
    proc, c(SLURM_MEM_PER_CPU = "3000", SLURM_CPUS_ON_NODE = "20"), native
  )
  expect_equal(unname(per_cpu$signals["slurm_allocation_remaining"]),
               58818560000)
  per_cpu_tasks <- mgcvST:::.mgcvst_memory_probe(
    proc, c(SLURM_MEM_PER_CPU = "3000", SLURM_CPUS_ON_NODE = "20",
            SLURM_NTASKS_PER_NODE = "1,4(x2)"), native
  )
  expect_equal(unname(per_cpu_tasks$signals["slurm_allocation_remaining"]),
               11632640000)
  ambiguous <- mgcvST:::.mgcvst_memory_probe(
    proc, c(SLURM_MEM_PER_CPU = "3000", SLURM_CPUS_ON_NODE = "20(x2)"), native
  )
  expect_false("slurm_allocation_remaining" %in% names(ambiguous$signals))
})

test_that("cgroup headroom wins over host and Slurm signals", {
  root <- tempfile("memory-priority-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  memory_fixture_write(root, "cg/job/memory.max", "60000000000")
  memory_fixture_write(root, "cg/job/memory.current", "20000000000")
  proc <- memory_proc_fixture(
    root, available_kb = 200000000, rss_kb = 1000000,
    cgroup = "0::/job",
    mountinfo = paste0("29 23 0:26 / ", file.path(root, "cg"),
                       " rw - cgroup2 cgroup rw")
  )
  z <- mgcvST:::.mgcvst_memory_probe(
    proc, c(SLURM_MEM_PER_NODE = "150000", SLURM_NTASKS_PER_NODE = "2"),
    c(available = NA_real_, resident = NA_real_)
  )
  expect_equal(z$available, 40000000000)
  expect_match(z$source, "cgroup2")
})

test_that("Windows-native readings are 64-bit inputs and unknown stays unknown", {
  root <- tempfile("memory-unknown-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  win <- mgcvST:::.mgcvst_memory_probe(
    root, character(),
    c(available = 120000000000, resident = 9000000000)
  )
  expect_equal(win$available, 120000000000)
  expect_equal(win$resident, 9000000000)
  unknown <- mgcvST:::.mgcvst_memory_probe(
    root, character(), c(available = NA_real_, resident = NA_real_)
  )
  expect_true(is.na(unknown$available))
  expect_true(is.na(unknown$resident))
  expect_identical(unknown$source, "unknown")
  expect_length(unknown$signals, 0L)
  zero <- mgcvST:::.mgcvst_memory_probe(
    root, character(), c(available = 0, resident = 0)
  )
  expect_identical(zero$available, 0)
  expect_identical(zero$resident, 0)
})
