# qc_utils.R
# Shared helpers for all qc_*() functions: one standard format for QC flags.

# Empty flag table (used when a check finds nothing)
empty_flags <- function() {
  data.frame(stage = character(), check = character(), severity = character(),
             subject = character(), session = character(), region = character(),
             value = numeric(), message = character())
}

# Make flags for the rows of `df` where `hit` is TRUE
make_flags <- function(df, hit, stage, check, severity, value, message) {
  hit[is.na(hit)] <- FALSE
  if (!any(hit)) return(empty_flags())
  data.frame(
    stage    = stage,
    check    = check,
    severity = severity,
    subject  = df$subject[hit],
    session  = df$session[hit],
    region   = if ("region" %in% names(df)) df$region[hit] else NA_character_,
    value    = as.numeric(value[hit]),
    message  = message
  )
}

# Add qc_exclude / qc_reason columns to the data from a flag table
apply_flags <- function(df, flags, keys) {
  df$qc_exclude <- FALSE
  df$qc_reason  <- ""
  ex <- flags[flags$severity == "exclude", ]
  for (i in seq_len(nrow(ex))) {
    row <- Reduce(`&`, lapply(keys, function(k) df[[k]] == ex[[k]][i]))
    df$qc_exclude[row] <- TRUE
    df$qc_reason[row]  <- trimws(paste(df$qc_reason[row], ex$check[i]))
  }
  df
}