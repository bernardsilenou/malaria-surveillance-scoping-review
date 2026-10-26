# Review ASySD's possible duplicate pairs, merge confirmed matches, and export
# the final records for title/abstract screening.
# Run this only after 01_deduplicate_and_audit_search_records.R.
# ASySD's documented manual review opens its pair-review Shiny interface.

if (!requireNamespace("ASySD", quietly = TRUE)) {
  stop("ASySD is not installed in this R library.")
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_dir, "databases", "deduplication_asysd")
result_path <- file.path(output_dir, "asysd_deduplication_result.rds")
if (!file.exists(result_path)) {
  stop("First run 01_deduplicate_and_audit_search_records.R to create the ASySD results file.")
}

if (as.character(utils::packageVersion("ASySD")) != "0.4.7") {
  stop("This review step was checked against ASySD 0.4.7. Review it before using a different version.")
}

saved <- readRDS(result_path)
dedup_result <- saved$dedup_result
source_records <- saved$source_records
if (is.null(dedup_result$unique) || is.null(source_records)) {
  stop("The saved ASySD result does not include the expected unique records and source records.")
}

manual_pairs <- dedup_result$manual_dedup
if (is.null(manual_pairs)) manual_pairs <- data.frame()
review_checkpoint_path <- file.path(output_dir, "manually_reviewed_possible_duplicate_pairs.csv")

if (nrow(manual_pairs) > 0) {
  # Reuse any earlier decisions so an interrupted or intentionally partial
  # review can be resumed. ASySD displays only rows whose result is still NA.
  review_input <- manual_pairs
  if (file.exists(review_checkpoint_path)) {
    saved_review <- utils::read.csv(
      review_checkpoint_path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA")
    )
    has_pair_ids <- all(c("record_id1", "record_id2") %in% names(saved_review))
    if (has_pair_ids) {
      pair_key <- function(data) {
        first_id <- as.character(data$record_id1)
        second_id <- as.character(data$record_id2)
        paste(pmin(first_id, second_id), pmax(first_id, second_id), sep = "::")
      }
      current_keys <- pair_key(manual_pairs)
      saved_keys <- pair_key(saved_review)
      if (!anyDuplicated(saved_keys)) {
        saved_rows <- match(current_keys, saved_keys)
        reused <- !is.na(saved_rows)
        if (any(reused)) {
          review_input <- manual_pairs
          review_input$result <- rep(NA_character_, nrow(manual_pairs))
          if ("result" %in% names(saved_review)) {
            review_input$result[reused] <- as.character(saved_review$result[saved_rows[reused]])
          }
          message(
            "Restored decisions for ", sum(reused), " of ", nrow(manual_pairs),
            " current pairs; ", sum(!reused), " new pair(s) need review."
          )
        } else {
          warning("No saved decisions match the current candidate pairs; starting a fresh review.")
        }
      } else {
        warning("Ignoring the saved manual-review file because it contains duplicate pair IDs.")
      }
    } else {
      warning("Ignoring the saved manual-review file because it lacks the expected pair IDs.")
    }
  }

  # ASySD marks each candidate as match, no_match, or NA (skip/unresolved).
  if (!"result" %in% names(review_input) || anyNA(review_input$result)) {
    reviewed_pairs <- ASySD::manual_dedup_shiny(review_input)
  } else {
    reviewed_pairs <- review_input
    message("All candidate pairs already have decisions; reusing the saved review.")
  }
  utils::write.csv(
    reviewed_pairs,
    review_checkpoint_path,
    row.names = FALSE,
    na = ""
  )

  if (!"result" %in% names(reviewed_pairs)) {
    stop("The ASySD review did not return a result column; no final records were created.")
  }
  confirmed_pairs <- dplyr::filter(reviewed_pairs, result == "match")
  if (nrow(confirmed_pairs) > 0) {
    final_unique <- ASySD::dedup_citations_add_manual(
      unique_citations = dedup_result$unique,
      merge_citations = TRUE,
      additional_pairs = confirmed_pairs,
      extra_merge_fields = c("source_record_id", "original_record_id", "accession_number", "import_file", "keywords")
    )
  } else {
    # ASySD 0.4.7 renames duplicate_id before returning early when no pairs
    # are marked as matches. Keep the valid automatic result unchanged here.
    final_unique <- dedup_result$unique
  }
} else {
  reviewed_pairs <- data.frame()
  final_unique <- dedup_result$unique
  utils::write.csv(
    reviewed_pairs,
    file.path(output_dir, "manually_reviewed_possible_duplicate_pairs.csv"),
    row.names = FALSE
  )
}

# Build a record-level audit from the final duplicate groups.
group_members <- lapply(seq_len(nrow(final_unique)), function(i) {
  ids <- final_unique$record_ids[[i]]
  if (is.na(ids) || !nzchar(ids)) ids <- final_unique$duplicate_id[[i]]
  strsplit(ids, ",\\s*")[[1]]
})
group_map <- dplyr::bind_rows(lapply(seq_along(group_members), function(i) {
  ids <- group_members[[i]]
  data.frame(
    record_id = ids,
    duplicate_group_id = final_unique$duplicate_id[[i]],
    dedup_status = ifelse(ids == final_unique$duplicate_id[[i]], "retained_representative", "duplicate_removed_after_review"),
    stringsAsFactors = FALSE
  )
}))
final_audit <- dplyr::left_join(source_records, group_map, by = "record_id")
if (any(is.na(final_audit$dedup_status))) {
  stop("Some imported records are missing from the final ASySD group audit.")
}
auto_representatives <- as.character(dedup_result$unique$duplicate_id)
final_audit$dedup_status <- ifelse(
  final_audit$dedup_status == "retained_representative",
  "retained_representative",
  ifelse(
    final_audit$record_id %in% auto_representatives,
    "duplicate_removed_after_manual_review",
    "automatically_removed_duplicate"
  )
)

utils::write.csv(
  final_audit,
  file.path(output_dir, "final_duplicate_group_audit_all_records.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  dplyr::filter(final_audit, dedup_status != "retained_representative"),
  file.path(output_dir, "final_removed_duplicate_records.csv"),
  row.names = FALSE,
  na = ""
)

# CSV is suitable for ASReview/prismAId; RIS is convenient for reference
# managers and is also accepted by ASReview.
screening_csv <- final_unique |>
  dplyr::transmute(
    id = duplicate_id,
    title = title,
    abstract = abstract,
    authors = author,
    year = year,
    doi = doi,
    journal = journal,
    keywords = keywords,
    source = source
  )
utils::write.csv(
  screening_csv,
  file.path(output_dir, "final_unique_records_for_screening.csv"),
  row.names = FALSE,
  na = ""
)
ASySD::write_citations(
  final_unique,
  type = "ris",
  filename = file.path(output_dir, "final_unique_records_for_screening.ris")
)
saveRDS(final_unique, file.path(output_dir, "asysd_final_unique_records.rds"))

manual_matches <- if ("result" %in% names(reviewed_pairs)) sum(reviewed_pairs$result == "match", na.rm = TRUE) else 0L
unresolved_pairs <- if ("result" %in% names(reviewed_pairs)) {
  sum(is.na(reviewed_pairs$result) | !nzchar(reviewed_pairs$result))
} else {
  0L
}
automatic_duplicates_removed <- nrow(source_records) - nrow(dedup_result$unique)
manual_duplicates_removed <- nrow(dedup_result$unique) - nrow(final_unique)
total_duplicates_removed <- nrow(source_records) - nrow(final_unique)
summary_lines <- c(
  paste("ASySD version:", as.character(utils::packageVersion("ASySD"))),
  paste("Date finalized:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Records entering deduplication:", nrow(source_records)),
  paste("Records after automatic deduplication:", nrow(dedup_result$unique)),
  paste("Automatically removed duplicate records:", automatic_duplicates_removed),
  paste("Unique records after automated and manual deduplication:", nrow(final_unique)),
  paste("Additional duplicates removed after manual review:", manual_duplicates_removed),
  paste("Total duplicate records removed:", total_duplicates_removed),
  paste("Possible pairs manually confirmed as duplicates:", manual_matches),
  paste("Unresolved possible duplicate pairs:", unresolved_pairs),
  "The screening CSV/RIS contain deduplicated records; duplicate removal is not eligibility exclusion."
)
writeLines(summary_lines, file.path(output_dir, "final_run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
