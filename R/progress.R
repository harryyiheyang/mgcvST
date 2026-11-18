# Progress lines of the slow steps (estimation chunks and pair evaluation).
#
# A line is printed with message() each time another tenth of the work is done,
# and at least every `heartbeat` seconds whenever work completes, so that a run
# whose first tenth is slow still shows that it is alive:
#
#   [null] 4,480 / 11,184 features (40%) | 12m 05s elapsed | ETA 18m 10s | 6.2 features/s
#
# Work restored from a checkpoint counts as done but not towards the rate.
# Plain lines (no carriage returns) keep batch logs readable.
# options(mgcvST.progress = FALSE) silences them.

.mgcvst_progress_enabled <- function() {
  isTRUE(getOption("mgcvST.progress", TRUE))
}

.mgcvst_format_count <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

.mgcvst_format_seconds <- function(s) {
  if (!is.finite(s)) return("?")
  s <- max(0, round(s))
  if (s < 60) return(sprintf("%ds", s))
  if (s < 3600) return(sprintf("%dm %02ds", s %/% 60, s %% 60))
  sprintf("%dh %02dm", s %/% 3600, (s %% 3600) %/% 60)
}

# A progress tracker of `total` units of work, `resumed` of them already done.
# Returns a list of functions: add(units) after a piece completes (fresh =
# FALSE for a piece restored from a checkpoint, which does not count towards
# the rate), done() at the end (prints the final line once).
.mgcvst_progress <- function(label, total, unit = "features", resumed = 0,
                             heartbeat = 60, enabled = .mgcvst_progress_enabled()) {
  state <- new.env(parent = emptyenv())
  state$done <- resumed
  state$fresh <- 0
  state$t0 <- proc.time()[["elapsed"]]
  state$last <- state$t0
  state$tenth <- if (total > 0) floor(10 * resumed / total) else 10
  state$finished <- FALSE
  emit <- function() {
    now <- proc.time()[["elapsed"]]
    elapsed <- now - state$t0
    rate <- if (elapsed > 0) state$fresh / elapsed else NA_real_
    left <- total - state$done
    eta <- if (left <= 0) 0 else if (is.finite(rate) && rate > 0) left / rate else NA_real_
    message(sprintf(
      "[%s] %s / %s %s (%d%%) | %s elapsed | ETA %s | %s %s/s",
      label, .mgcvst_format_count(state$done), .mgcvst_format_count(total),
      unit, if (total > 0) as.integer(floor(100 * state$done / total)) else 100L,
      .mgcvst_format_seconds(elapsed), .mgcvst_format_seconds(eta),
      if (is.finite(rate)) format(signif(rate, 3), big.mark = ",") else "?", unit
    ))
    state$last <- now
  }
  if (enabled && total > 0) {
    message(sprintf("[%s] %s %s%s", label, .mgcvst_format_count(total), unit,
                    if (resumed > 0) sprintf(", %s restored from the checkpoint",
                                             .mgcvst_format_count(resumed)) else ""))
  }
  add <- function(units, fresh = TRUE) {
    if (length(units) != 1L || !is.finite(units)) return(invisible(NULL))
    state$done <- min(total, state$done + units)
    if (fresh) state$fresh <- state$fresh + units
    if (!enabled || total <= 0) return(invisible(NULL))
    tenth <- floor(10 * state$done / total)
    now <- proc.time()[["elapsed"]]
    if (tenth > state$tenth || now - state$last >= heartbeat) {
      state$tenth <- max(state$tenth, tenth)
      if (state$done >= total) state$finished <- TRUE
      emit()
    }
    invisible(NULL)
  }
  done <- function() {
    if (enabled && total > 0 && !state$finished) {
      state$finished <- TRUE
      emit()
    }
    invisible(NULL)
  }
  list(add = add, done = done)
}
