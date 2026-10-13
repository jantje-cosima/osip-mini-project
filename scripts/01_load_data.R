# 01_load_data.R
# Stage 1 (Load): read input files and check their STRUCTURE
# (columns, types, keys). Plausibility of the values is checked
# later in 02_qc_input.R.

# --- Helpers -----------------------------------------------------------------

# Stop if any required column is missing
check_columns <- function(df, required, file) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(file, ": missing column(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# Stop if any of the given columns fails a type check (e.g. is.numeric)
check_types <- function(df, cols, type_check, type_name, file) {
  bad <- cols[!vapply(df[cols], type_check, logical(1))]
  if (length(bad) > 0) {
    stop(file, ": column(s) not ", type_name, ": ", paste(bad, collapse = ", "))
  }
  invisible(TRUE)
}

# Stop if the key columns do not identify each row uniquely
check_unique_keys <- function(df, keys, file) {
  dup <- duplicated(df[keys])
  if (any(dup)) {
    stop(file, ": ", sum(dup), " duplicated row(s) for key ",
         paste(keys, collapse = " x "))
  }
  invisible(TRUE)
}

# --- VT estimates --------------------------------------------------------------

load_vt <- function(path = file.path("data", "raw", "vt_estimates.csv")) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.csv(path, stringsAsFactors = FALSE)

  keys    <- c("subject", "session", "region")
  numeric <- c("K1", "k2", "k3", "k4", "VT", "VT_SE", "VT_SE_pct", "wrss")
  logical <- c("converged", "at_bound")

  check_columns(df, c(keys, numeric, logical), path)
  check_types(df, numeric, is.numeric, "numeric", path)
  check_types(df, logical, is.logical, "logical (TRUE/FALSE)", path)
  check_unique_keys(df, keys, path)

  bad_session <- setdiff(unique(df$session), c("baseline", "block1", "block2"))
  if (length(bad_session) > 0) {
    stop(path, ": unknown session(s): ", paste(bad_session, collapse = ", "))
  }

  df
}

# --- Plasma drug concentrations ------------------------------------------------

load_drug <- function(path = file.path("data", "raw", "drug_plasma.csv")) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.csv(path, stringsAsFactors = FALSE)

  keys    <- c("subject", "session", "sample")
  numeric <- c("dose_mg", "time_post_dose_h", "drug_ng_mL")

  check_columns(df, c(keys, numeric), path)
  check_types(df, numeric, is.numeric, "numeric", path)
  check_unique_keys(df, keys, path)

  bad_session <- setdiff(unique(df$session), c("block1", "block2"))
  if (length(bad_session) > 0) {
    stop(path, ": unexpected session(s): ", paste(bad_session, collapse = ", "))
  }

  df
}