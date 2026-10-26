# Deduplicate database search records with ASySD and export screening files.
# Run this script from the sormas-malaria R project root.
# The source exports in databases/ are read-only inputs; all outputs go to
# databases/deduplication_asysd/.
#
# ASySD 0.4.7's non-Shiny wrapper omits the required extra_merge_fields
# argument when it calls order_citations(). The compatibility wrapper below
# follows ASySD's non-Shiny workflow and supplies that argument. It calls
# ASySD's matching, classification, and merge functions; no package files are
# modified. Review this wrapper if changing the installed ASySD version.

if (!requireNamespace("ASySD", quietly = TRUE)) {
  stop("ASySD is not installed in this R library. Install it before running this script.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("The dplyr package is required; it is normally installed with ASySD.")
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_dir, "databases")
output_dir <- file.path(data_dir, "deduplication_asysd")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_started <- Sys.time()
run_id <- format(run_started, "%Y%m%d_%H%M%S")
log_path <- file.path(output_dir, paste0("deduplication_run_", run_id, ".log"))
log_progress <- function(..., .newline = TRUE) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", paste0(..., collapse = ""))
  cat(line, if (.newline) "\n" else "", sep = "")
  cat(line, if (.newline) "\n" else "", sep = "", file = log_path, append = TRUE)
  flush.console()
  invisible(line)
}

log_progress("Run started. Working directory: ", project_dir)
log_progress("ASySD version: ", as.character(utils::packageVersion("ASySD")))

required_files <- c(
  "pubmed.txt",
  "Scopus.ris",
  "Web_of_Science_1.ris",
  "Web_of_Science_2.ris",
  "CINHAL.ris",
  "My_Library_Zotero_combined.ris"
)
missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0) {
  stop("Missing input file(s) in databases/: ", paste(missing_files, collapse = ", "))
}

# Convert PubMed's NBIB export to the common ASySD field structure.
read_pubmed_nbib <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  records <- list()
  current <- list()
  current_tag <- NULL

  add_record <- function() {
    if (length(current) == 0) return(invisible(NULL))
    get_values <- function(tag) unlist(current[[tag]], use.names = FALSE)
    get_one <- function(tag) {
      value <- get_values(tag)
      if (length(value) == 0) NA_character_ else paste(value, collapse = " ")
    }
    extract_doi <- function(value) {
      if (is.na(value) || !nzchar(value)) return(NA_character_)
      hits <- regmatches(value, gregexpr("10\\.[0-9]{4,9}/[^[:space:]]+", value, ignore.case = TRUE))[[1]]
      if (length(hits) == 0 || identical(hits, "")) return(NA_character_)
      tolower(sub("[.,;]+$", "", hits[[1]]))
    }

    pmid <- get_one("PMID")
    doi_text <- paste(c(get_values("LID"), get_values("AID")), collapse = " ")
    dp <- get_one("DP")
    year_match <- if (is.na(dp)) character(0) else regmatches(dp, regexpr("(19|20)[0-9]{2}", dp))
    year <- if (length(year_match) == 0) NA_character_ else year_match[[1]]
    next_record_number <- length(records) + 1L

    title <- get_one("TI")
    if (is.na(title) || !nzchar(title)) title <- get_one("BTI")

    records[[next_record_number]] <<- data.frame(
      title = title,
      author = if (length(get_values("AU")) == 0) NA_character_ else paste(get_values("AU"), collapse = "; "),
      year = as.character(year),
      journal = if (!is.na(get_one("JT"))) get_one("JT") else get_one("TA"),
      doi = extract_doi(doi_text),
      pages = get_one("PG"),
      volume = get_one("VI"),
      number = get_one("IP"),
      abstract = get_one("AB"),
      keywords = if (length(get_values("OT")) == 0) NA_character_ else paste(get_values("OT"), collapse = "; "),
      record_id = paste0("PM_", sprintf("%05d", next_record_number)),
      source_record_id = pmid,
      original_record_id = pmid,
      accession_number = pmid,
      label = NA_character_,
      source = "PubMed",
      import_file = basename(path),
      url = NA_character_,
      isbn = get_one("IS"),
      stringsAsFactors = FALSE
    )
    current <<- list()
    current_tag <<- NULL
    invisible(NULL)
  }

  for (line in lines) {
    field_match <- regmatches(line, regexec("^([A-Z][A-Z0-9]{1,3})\\s*-\\s*(.*)$", line))[[1]]
    if (length(field_match) > 0) {
      tag <- field_match[[2]]
      value <- field_match[[3]]
      if (tag == "PMID" && length(current) > 0) add_record()
      current_tag <- tag
      current[[tag]] <- c(current[[tag]], value)
    } else if (grepl("^\\s+\\S", line) && !is.null(current_tag)) {
      values <- current[[current_tag]]
      values[[length(values)]] <- paste(values[[length(values)]], trimws(line))
      current[[current_tag]] <- values
    }
  }
  add_record()

  if (length(records) == 0) stop("No PubMed records were parsed from ", path)
  dplyr::bind_rows(records)
}

read_ris_source <- function(file_name, source_name, id_prefix) {
  path <- file.path(data_dir, file_name)
  records <- ASySD::load_search(path, method = "ris")
  if (nrow(records) == 0) stop("No records were parsed from ", file_name)

  original_id <- if ("record_id" %in% names(records)) as.character(records$record_id) else rep(NA_character_, nrow(records))
  accession <- if ("AN" %in% names(records)) as.character(records$AN) else rep(NA_character_, nrow(records))
  source_id <- original_id
  missing_id <- is.na(source_id) | !nzchar(source_id)
  source_id[missing_id] <- accession[missing_id]

  records$original_record_id <- original_id
  records$source_record_id <- source_id
  records$accession_number <- accession
  records$keywords <- if ("KW" %in% names(records)) as.character(records$KW) else NA_character_
  records$record_id <- paste0(id_prefix, "_", sprintf("%05d", seq_len(nrow(records))))
  records$source <- source_name
  records$import_file <- file_name
  if (!"label" %in% names(records)) records$label <- NA_character_
  records
}

log_progress("Importing source exports.")
pubmed <- read_pubmed_nbib(file.path(data_dir, "pubmed.txt"))
scopus <- read_ris_source("Scopus.ris", "Scopus", "SC")
wos_1 <- read_ris_source("Web_of_Science_1.ris", "Web of Science", "W1")
wos_2 <- read_ris_source("Web_of_Science_2.ris", "Web of Science", "W2")
cinahl <- read_ris_source("CINHAL.ris", "CINAHL", "CI")

source_records <- dplyr::bind_rows(pubmed, scopus, wos_1, wos_2, cinahl)
source_records$record_id <- as.character(source_records$record_id)
if (anyDuplicated(source_records$record_id)) stop("Generated record IDs are not unique.")
log_progress("Imported ", nrow(source_records), " records from ", length(unique(source_records$source)), " source databases.")

# Standardize columns expected by ASySD and retain provenance fields for audit.
required_columns <- c(
  "author", "year", "journal", "doi", "title", "pages", "volume", "number",
  "abstract", "record_id", "isbn", "label", "source", "source_record_id",
  "original_record_id", "accession_number", "import_file", "url", "keywords"
)
for (column in setdiff(required_columns, names(source_records))) source_records[[column]] <- NA_character_
source_records <- source_records[, unique(c(required_columns, names(source_records)[!names(source_records) %in% required_columns]))]

source_counts <- source_records |>
  dplyr::count(source, import_file, name = "records_imported") |>
  dplyr::arrange(source, import_file)
utils::write.csv(source_counts, file.path(output_dir, "source_counts.csv"), row.names = FALSE, na = "")
utils::write.csv(source_records, file.path(output_dir, "all_imported_records_with_source.csv"), row.names = FALSE, na = "")
saveRDS(source_records, file.path(output_dir, "all_imported_records_with_source.rds"))
log_progress("Saved import checkpoint: all_imported_records_with_source.rds")

# The combined Zotero export is used only for an import-count check; it is not
# added to the database exports, which would duplicate the search records.
zotero_combined <- ASySD::load_search(file.path(data_dir, "My_Library_Zotero_combined.ris"), method = "ris")
count_check <- data.frame(
  item = c("records in the four original database exports", "records in the combined Zotero RIS"),
  n = c(nrow(source_records), nrow(zotero_combined)),
  stringsAsFactors = FALSE
)
utils::write.csv(count_check, file.path(output_dir, "input_count_check.csv"), row.names = FALSE)
log_progress("Wrote input count check. Original source exports: ", nrow(source_records), "; combined Zotero RIS: ", nrow(zotero_combined), ".")

export_results <- function(dedup_result) {
  unique_records <- dedup_result$unique
  if (is.null(unique_records) || nrow(unique_records) == 0) {
    stop("ASySD did not return a unique-record dataset.")
  }

  # Export a tabular file for ASReview/prismAId and RIS for reference managers.
  asreview_csv <- unique_records |>
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
  utils::write.csv(asreview_csv, file.path(output_dir, "auto_unique_pre_manual_review.csv"), row.names = FALSE, na = "")
  ASySD::write_citations(unique_records, type = "ris", filename = file.path(output_dir, "auto_unique_pre_manual_review.ris"))

  # Expand ASySD groups to one row per original record for the audit trail.
  group_members <- lapply(seq_len(nrow(unique_records)), function(i) {
    ids <- unique_records$record_ids[[i]]
    if (is.na(ids) || !nzchar(ids)) ids <- unique_records$duplicate_id[[i]]
    strsplit(ids, ",\\s*")[[1]]
  })
  group_map <- dplyr::bind_rows(lapply(seq_along(group_members), function(i) {
    ids <- group_members[[i]]
    data.frame(
      record_id = ids,
      duplicate_group_id = unique_records$duplicate_id[[i]],
      dedup_status = ifelse(ids == unique_records$duplicate_id[[i]], "retained_representative", "automatically_removed_duplicate"),
      stringsAsFactors = FALSE
    )
  }))
  duplicate_audit <- dplyr::left_join(source_records, group_map, by = "record_id")
  if (any(is.na(duplicate_audit$dedup_status))) {
    stop("Some imported records are missing from the ASySD group audit.")
  }

  utils::write.csv(duplicate_audit, file.path(output_dir, "duplicate_group_audit_all_records.csv"), row.names = FALSE, na = "")
  utils::write.csv(
    dplyr::filter(duplicate_audit, dedup_status == "automatically_removed_duplicate"),
    file.path(output_dir, "automatically_removed_duplicate_records.csv"),
    row.names = FALSE, na = ""
  )

  manual_pairs <- dedup_result$manual_dedup
  if (is.null(manual_pairs)) manual_pairs <- data.frame()
  utils::write.csv(manual_pairs, file.path(output_dir, "possible_duplicate_pairs_for_manual_review.csv"), row.names = FALSE, na = "")

  saveRDS(
    list(dedup_result = dedup_result, source_records = source_records),
    file.path(output_dir, "asysd_deduplication_result.rds")
  )
  summary_lines <- c(
    paste("ASySD version:", as.character(utils::packageVersion("ASySD"))),
    paste("Date run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Records imported from the four original database exports:", nrow(source_records)),
    paste("Records in the combined Zotero RIS (count check only):", nrow(zotero_combined)),
    paste("Provisional records after automatic deduplication:", nrow(unique_records)),
    paste("Automatically removed duplicate records:", sum(duplicate_audit$dedup_status == "automatically_removed_duplicate")),
    paste("Possible duplicate pairs requiring manual review:", nrow(manual_pairs)),
    "Automatic duplicate removal is record deduplication, not title/abstract eligibility exclusion.",
    "The CSV/RIS outputs are provisional until possible_duplicate_pairs_for_manual_review.csv has been reviewed with ASySD's manual review step.",
    "Use 02_review_possible_duplicate_pairs.R to review candidate pairs and export the final screening set."
  )
  writeLines(summary_lines, file.path(output_dir, "run_summary.txt"))
  summary_lines
}

asysd_dedup_console <- function(
    raw_citations,
    manual_dedup = TRUE,
    merge_citations = TRUE,
    keep_source = NULL,
    keep_label = NULL,
    extra_merge_fields = NULL,
    show_unknown_tags = TRUE,
    user_input = 1) {

  # These are ASySD internal functions. The version check makes the dependency
  # on the verified package implementation explicit and prevents silent drift.
  if (as.character(utils::packageVersion("ASySD")) != "0.4.7") {
    stop("This compatibility wrapper was checked against ASySD 0.4.7. Review it before using a different version.")
  }

  required_internal <- c(
    "format_rerun", "add_missing_cols", "order_citations", "format_citations",
    "match_citations", "identify_true_matches", "generate_dup_id",
    "merge_metadata", "keep_one_unique_citation", "process_possible_pairs"
  )
  missing_internal <- required_internal[!vapply(required_internal, exists, logical(1),
    envir = asNamespace("ASySD"), inherits = FALSE
  )]
  if (length(missing_internal) > 0) {
    stop("ASySD internal function(s) not found: ", paste(missing_internal, collapse = ", "))
  }

  raw_citations <- ASySD:::format_rerun(raw_citations)
  cols <- c("author", "year", "journal", "doi", "title", "pages", "volume", "number", "abstract", "record_id", "isbn", "label", "source")
  essential_cols <- c("author", "year", "journal", "doi", "title", "pages", "volume", "number", "abstract", "record_id", "isbn")
  all_missing_cols <- cols[!(cols %in% names(raw_citations))]
  missing_cols <- essential_cols[!(essential_cols %in% names(raw_citations))]

  if (length(missing_cols) > 0 && is.na(user_input)) {
    user_input <- utils::menu(c("Yes", "No"), title = "Proceed with missing citation fields?")
  }
  if (length(missing_cols) > 0 && user_input != 1) return("Halting dedup...")
  if (length(missing_cols) > 0) {
    message("Proceeding with missing fields: ", paste(missing_cols, collapse = ", "))
  }
  raw_citations[all_missing_cols] <- NA

  if (show_unknown_tags) {
    raw_citations <- raw_citations |>
      dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(., ""))) |>
      dplyr::mutate(
        label = ifelse(is.na(.data$label), "unknown", .data$label),
        source = ifelse(is.na(.data$source), "unknown", .data$source)
      )
    if (length(extra_merge_fields) > 0) {
      raw_citations <- raw_citations |>
        dplyr::mutate(dplyr::across(dplyr::all_of(extra_merge_fields), ~ tidyr::replace_na(., "unknown")))
    }
  } else {
    raw_citations <- raw_citations |>
      dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(., "")))
  }

  raw_citations <- ASySD:::add_missing_cols(raw_citations)
  # Custom provenance fields are preserved in raw_citations and are merged
  # after matching. They are not inputs to ASySD's similarity scores, so do
  # not pass them through the older order_citations() selector.
  ordered_citations <- ASySD:::order_citations(raw_citations, NULL)
  formatted_citations <- ASySD:::format_citations(ordered_citations)
  pairs <- ASySD:::match_citations(formatted_citations)

  result <- if (manual_dedup) list(manual_dedup = data.frame()) else list()
  if (is.null(pairs) || length(pairs$record_id1) == 0) {
    warning("No duplicate candidates detected.")
    raw_citations$duplicate_id <- raw_citations$record_ids <- raw_citations$record_id
    result$unique <- raw_citations
    return(result)
  }

  pair_types <- ASySD:::identify_true_matches(pairs)
  true_pairs <- pair_types$true_pairs
  if (is.null(true_pairs) || length(true_pairs$record_id1) == 0) {
    warning("No automatically confirmed duplicate citations detected.")
    raw_citations$duplicate_id <- raw_citations$record_ids <- raw_citations$record_id
    result$unique <- raw_citations
    if (manual_dedup) {
      no_auto_ids <- raw_citations |>
        dplyr::transmute(record_id = as.character(.data$record_id), duplicate_id = as.character(.data$record_id))
      result$manual_dedup <- ASySD:::process_possible_pairs(
        pair_types$maybe_pairs, ordered_citations, no_auto_ids,
        # ASySD 0.4.7 builds these dynamic column names with `:=` in a way
        # that fails for character vectors under current rlang/dplyr. These
        # extra fields are for displaying candidate pairs; standard fields,
        # including source and record IDs, are sufficient for review. The
        # fields remain in source_records and are merged into final records.
        NULL
      )
    }
    return(result)
  }

  matched_pairs <- ASySD:::generate_dup_id(true_pairs, raw_citations, keep_source, keep_label)
  if (merge_citations) {
    result$unique <- ASySD:::merge_metadata(matched_pairs, extra_merge_fields)
  } else {
    result$unique <- ASySD:::keep_one_unique_citation(matched_pairs)
  }
  if (manual_dedup) {
    result$manual_dedup <- ASySD:::process_possible_pairs(
      pair_types$maybe_pairs, ordered_citations, matched_pairs,
      # Avoid ASySD 0.4.7's incompatible dynamic `:=` handling for custom
      # display columns. Source and record IDs are still shown in pair review.
      NULL
    )
  }
  result$unique <- dplyr::ungroup(result$unique)

  n_start <- length(unique(formatted_citations$record_id))
  n_unique <- length(unique(result$unique$duplicate_id))
  message(n_start, " citations loaded; ", n_start - n_unique,
    " duplicate citations removed; ", n_unique, " unique citations remaining.")
  result
}

log_progress("Starting ASySD automatic matching. This can take time; the import checkpoint is already saved.")
dedup_result <- tryCatch(
  asysd_dedup_console(
  raw_citations = source_records,
  manual_dedup = TRUE,
  merge_citations = TRUE,
  extra_merge_fields = c("source_record_id", "original_record_id", "accession_number", "import_file", "keywords"),
  show_unknown_tags = TRUE,
  user_input = 1
  ),
  error = function(e) {
    log_progress("ASySD automatic matching failed: ", conditionMessage(e))
    stop(e)
  }
)
log_progress("ASySD automatic matching completed. Unique groups: ", nrow(dedup_result$unique), "; possible pairs for manual review: ", if (is.null(dedup_result$manual_dedup)) 0L else nrow(dedup_result$manual_dedup), ".")
summary_lines <- export_results(dedup_result)
log_progress("Automatic-stage outputs and ASySD result checkpoint saved successfully.")
message(paste(summary_lines, collapse = "\n"))
