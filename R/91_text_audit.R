## =====================================================================
## 91_text_audit.R -- do the verified numbers actually appear in the text?
## =====================================================================
## The verification scripts check that each number a paper reports can be
## recomputed from the stored output. They do not check the other
## direction: that the manuscript still contains the value that was
## verified. A number updated in the pipeline but left stale in the text
## passes verification and reaches the reviewer.
##
## This script closes that gap for both papers. It reads the claim table
## each verification script writes, and for every claim marked as passing
## it searches the manuscript for the value, formatted as a paper would
## format it. Values it cannot find are listed: either the text is stale,
## or it rounds them differently, or the claim is not quoted in the text
## at all (a table entry checked but not repeated in prose, say). Each
## case needs a human decision, so the script reports rather than judges.
##
## It also counts \TBD{} markers, which must be zero before submission.
##
## Run from the repository root:  source("R/91_text_audit.R")
## 90_verify_claims_B.R already contains this layer for Paper B; running
## it here as well is a harmless cross-check.
## =====================================================================

## paper -> claim table, candidate manuscript paths
## Set TEX_A / TEX_B before sourcing to name the manuscripts explicitly.
## Do that if several copies exist: the search below takes the first match,
## which may be an old one.
JOBS <- list(
  A = list(csv = "results/claim_verification.csv",
           tex = c(if (exists("TEX_A")) TEX_A, "PaperA.tex", "../PaperA.tex")),
  B = list(csv = "results/real/verify_claims_B.csv",
           tex = c(if (exists("TEX_B")) TEX_B, "PaperB.tex", "../PaperB.tex"))
)

find_tex <- function(cand, pattern) {
  hit <- cand[file.exists(cand)]
  if (!length(hit))
    hit <- list.files(".", pattern, recursive = TRUE, full.names = TRUE)
  if (!length(hit))
    hit <- list.files("..", pattern, recursive = TRUE, full.names = TRUE)
  if (length(hit)) hit[1] else NA_character_
}

## which column holds the value as the paper states it, and which the status
pick <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) hit[1] else NA_character_
}

## the ways a number of this size is plausibly written in a manuscript
## A paper rounds to as few digits as the quantity deserves, so we try
## several roundings -- but only those that still identify the number.
## Rounding 0.087 to one digit gives 0.1, which matches almost any text,
## so coarse roundings are allowed only for values large enough to remain
## distinctive: whole numbers from 10 up, one decimal from 1 up.
formats <- function(x) {
  if (!is.finite(x)) return(character(0))
  a <- abs(x)
  f <- c(sprintf("%.2f", x), sprintf("%.3f", x), sprintf("%.4f", x))
  if (a >= 1)  f <- c(f, sprintf("%.1f", x))
  if (a >= 10) f <- c(f, sprintf("%.0f", x))
  f <- c(f, sub("^(-?)0\\.", "\\1.", f))          # papers often drop the 0
  unique(f[nzchar(f)])
}

## match a number as a whole token: "1.17" must not be found inside
## "1.175" or "11.17"
in_text <- function(p, tex) {
  ## The only regex metacharacter a formatted number can contain is ".".
  ## With fixed = TRUE the replacement is literal, so it must be exactly
  ## "\\." here: writing "\\\\." puts a real backslash into the pattern,
  ## which then matches nothing and reports every decimal as missing.
  ## the lookbehind also excludes a sign, so that searching for 0.09 does
  ## not succeed on a "-0.09" in the text
  rx <- paste0("(?<![0-9.\\-\u2212])", gsub(".", "\\.", p, fixed = TRUE),
               "(?![0-9])")
  grepl(rx, tex, perl = TRUE)
}

audit_one <- function(tag, job) {
  cat("\n", strrep("=", 78), "\nPAPER ", tag, "\n", strrep("=", 78), "\n", sep = "")
  if (!file.exists(job$csv)) {
    cat("claim table not found at ", job$csv,
        "\n  run the verification script for this paper first.\n", sep = "")
    return(invisible(NULL))
  }
  tex_path <- find_tex(job$tex, sprintf("^Paper%s1?\\.tex$", tag))
  if (is.na(tex_path)) {
    cat("manuscript not found; set it by hand in JOBS at the top.\n")
    return(invisible(NULL))
  }
  cat("claims     : ", job$csv, "\nmanuscript : ", tex_path, "\n", sep = "")

  cl  <- read.csv(job$csv, stringsAsFactors = FALSE)
  tex <- paste(readLines(tex_path, warn = FALSE), collapse = " ")

  v_col <- pick(names(cl), c("stated", "value", "expected", "reported", "paper"))
  s_col <- pick(names(cl), c("status", "result", "pass", "ok", "verdict"))
  if (is.na(v_col)) {
    cat("cannot tell which column holds the stated value. Columns are:\n  ",
        paste(names(cl), collapse = ", "),
        "\nadd its name to the candidates in pick().\n", sep = "")
    return(invisible(NULL))
  }
  lab <- pick(names(cl), c("claim", "what", "description", "parameter"))
  sec <- pick(names(cl), c("section", "where", "location"))

  ok <- if (is.na(s_col)) rep(TRUE, nrow(cl)) else
    grepl("ok|pass|true", cl[[s_col]], ignore.case = TRUE)
  cat(sprintf("%d claims, %d of them passing and searchable\n", nrow(cl), sum(ok)))

  sub <- cl[ok, , drop = FALSE]
  val <- suppressWarnings(as.numeric(sub[[v_col]]))
  ## Values reported to fewer digits give more matches by chance, so a hit
  ## on a one-digit rounding is weak evidence; it is still better than
  ## flagging every rounded number in the paper.
  found <- vapply(seq_along(val), function(i) {
    f <- formats(val[i])
    if (!length(f)) return(NA)
    any(vapply(f, \(p) in_text(p, tex), TRUE))
  }, logical(1))

  bad <- sub[!is.na(found) & !found, , drop = FALSE]
  if (nrow(bad)) {
    cat("\nverified values NOT found in the manuscript (check each):\n")
    keep <- c(sec, lab, v_col); keep <- keep[!is.na(keep)]
    print(bad[, keep, drop = FALSE], row.names = FALSE)
  } else {
    cat("\nevery verified value appears in the manuscript.\n")
  }

  n_tbd <- length(gregexpr("\\\\TBD\\{", tex)[[1]])
  n_tbd <- if (regexpr("\\\\TBD\\{", tex) < 0) 0L else n_tbd - 1L
  cat(sprintf("TBD markers: %d%s\n", n_tbd,
              if (n_tbd > 0) "  <- must be 0 before submission" else ""))
  invisible(bad)
}

invisible(lapply(names(JOBS), \(k) audit_one(k, JOBS[[k]])))

cat("\nIf a claim belongs to the OTHER paper, the verification script for this\n")
cat("one is still checking the combined claim set and should be split.\n")
cat("\nNote: a value may legitimately be absent from the text -- a table entry\n")
cat("that prose does not repeat, or one the paper rounds differently. The\n")
cat("script cannot tell those from a stale number, so each hit needs a look.\n")
