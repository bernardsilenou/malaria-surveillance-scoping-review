# Malaria surveillance requirements scoping review

This R project supports record deduplication and preparation of the literature set for screening in a scoping review mapping requirements for malaria surveillance.

## Workflow

1. Place the original database exports in `databases/` using the filenames expected by `01_deduplicate_and_audit_search_records.R`.
2. Open `sormas-malaria.Rproj` in RStudio and run `01_deduplicate_and_audit_search_records.R` to import the exports, retain source provenance, run ASySD, and produce the provisional unique set and duplicate audit.
3. Run `02_review_possible_duplicate_pairs.R` to review ASySD's candidate pairs and produce the final screening files. The script can resume from its saved manual-review decisions.

## Main outputs

Outputs are written to `databases/deduplication_asysd/`, including:

- `final_unique_records_for_screening.csv` and `.ris`
- `final_duplicate_group_audit_all_records.csv`
- `final_removed_duplicate_records.csv`
- `final_run_summary.txt`

The current workflow used ASySD 0.4.7. Automatic deduplication and manual review identify duplicate records; they do not make eligibility decisions.

## Data handling

Search exports, citation-level outputs, and the working protocol are excluded from Git by default. They may contain database-supplied metadata, abstracts, and unpublished review material. Keep the original exports backed up separately and only add these files to a repository if sharing is permitted and the repository access level is appropriate.
