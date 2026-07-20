#!/usr/bin/env Rscript
# =============================================================================
# build_ion_library.R
#
# Combines:
#   - existing ion_library.csv (D:/__STUDIES/ion_library.csv)
#   - all sheets in _method_panels.xlsx
# and joins metadata from:
#   - __featureMetadata_all.xlsx  (NEW oxylipin_metadata, SL_metadata,
#                                   PolarNeg_metadata, PolarPos_metadata)
#
# Output:
#   D:/__STUDIES/ion_library_combined.csv
#
# Output schema:
#   core:        transition_id, precursor_mz, product_mz, collision_eV,
#                cone_V, Polarity, Type, Lipid, Class
#   oxylipin:    PUFA, Primary, Enzymatic_pathway, Immediate_Pathway,
#                Autoxidation, Enzyme_1, Precursors, Omega_6, Omega_3
#   SL:          Class_abbrev, DeNovo_synthesis, Enzyme_list, GBA-GALC,
#                PPAP2A, ENPP7-SMPD, DEGS, CERK, CERS1..CERS6, UGT8-UGCG,
#                GLA-NEU, GLB-NEU, B4GALT, SGMS, SGMS1, SPHK, KDSR, SPTLC,
#                SGPP, ACER, Unkown
#   polar:       Super_Pathway, Sub_Pathway
#   provenance:  source (which panels), flag (ok / conflict / blank)
# =============================================================================

suppressMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# ---- inputs ----------------------------------------------------------------
LIB_CSV  <- "D:/__STUDIES/ion_library.csv"
PANELS   <- "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Panels/_method_panels.xlsx"
META     <- "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Panels/__featureMetadata_all.xlsx"
OUT_CSV  <- "D:/__STUDIES/ion_library_combined.csv"

# LC-MSMS feature_metadata sheets — lowercase `class` column is preferred
# source for the Class column when a transition matches by Processing_name.
LCMS_FILES <- c(
  "D:/__STUDIES/021_HDM_mice/data/LC-MSMS/PL_neg_mouseSections_HDM.xlsx",
  "D:/__STUDIES/021_HDM_mice/data/LC-MSMS/PL_pos_mouseSections_HDM.xlsx",
  "D:/__STUDIES/021_HDM_mice/data/LC-MSMS/polarPos_mouseSections_HDM.xlsx",
  "D:/__STUDIES/021_HDM_mice/data/LC-MSMS/SL_pos_mouseSections_HDM.xlsx"
)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 1. read existing ion library ------------------------------------------
base_lib <- read_csv(LIB_CSV, show_col_types = FALSE) |>
  rename(any_of(c(transition_id = "transition_id",
                  precursor_mz  = "precursor_mz",
                  product_mz    = "product_mz",
                  collision_eV  = "collision_eV",
                  cone_V        = "cone_V",
                  Polarity      = "Polarity",
                  Type          = "Type",
                  Lipid         = "Lipid"))) |>
  mutate(source = "ion_library.csv")

cat(sprintf("[base library]  %d rows\n", nrow(base_lib)))

# ---- 2. helper: pull a column safely (returns NA col if missing) -----------
col_or_na <- function(df, name, n = nrow(df)) {
  if (name %in% names(df)) df[[name]] else rep(NA, n)
}

# ---- 3. normalize one panel sheet to canonical schema ----------------------
normalize_panel <- function(sheet) {
  df <- suppressWarnings(suppressMessages(read_excel(PANELS, sheet)))
  if (nrow(df) == 0) return(NULL)
  nm <- names(df)

  pick_lipid <- function(d) {
    if      ("Lipid_class" %in% names(d))   d$Lipid_class
    else if ("Lipid class" %in% names(d))   d$`Lipid class`
    else                                    rep(NA_character_, nrow(d))
  }

  if ("transition_name" %in% nm) {
    # ---- Schema A (single MRM per row) ----
    out <- tibble(
      transition_id = as.character(df$transition_name),
      precursor_mz  = suppressWarnings(as.numeric(col_or_na(df, "Parent"))),
      product_mz    = suppressWarnings(as.numeric(col_or_na(df, "Daughter"))),
      collision_eV  = suppressWarnings(as.numeric(col_or_na(df, "CE"))),
      cone_V        = suppressWarnings(as.numeric(col_or_na(df, "Cone V"))),
      Polarity      = as.character(col_or_na(df, "Polarity")),
      Type          = as.character(col_or_na(df, "Type")),
      Lipid         = as.character(pick_lipid(df)),
      source        = sheet
    )
  } else if ("Process_name" %in% nm) {
    # ---- Schema B (1 or 2 MRMs per row) ----
    rows1 <- tibble(
      transition_id = as.character(df$Process_name),
      precursor_mz  = suppressWarnings(as.numeric(col_or_na(df, "Precursor1"))),
      product_mz    = suppressWarnings(as.numeric(col_or_na(df, "Product1"))),
      collision_eV  = suppressWarnings(as.numeric(col_or_na(df, "CE1"))),
      cone_V        = suppressWarnings(as.numeric(col_or_na(df, "Cone1"))),
      Polarity      = as.character(col_or_na(df, "Polarity")),
      Type          = as.character(col_or_na(df, "Type")),
      Lipid         = as.character(pick_lipid(df)),
      source        = sheet
    )
    rows2 <- if ("Process_name2" %in% nm) {
      tibble(
        transition_id = as.character(df$Process_name2),
        precursor_mz  = suppressWarnings(as.numeric(col_or_na(df, "Precursor2"))),
        product_mz    = suppressWarnings(as.numeric(col_or_na(df, "Product2"))),
        collision_eV  = suppressWarnings(as.numeric(col_or_na(df, "CE2"))),
        cone_V        = suppressWarnings(as.numeric(col_or_na(df, "Cone2"))),
        Polarity      = as.character(col_or_na(df, "Polarity")),
        Type          = as.character(col_or_na(df, "Type")),
        Lipid         = as.character(pick_lipid(df)),
        source        = sheet
      ) |> filter(!is.na(transition_id))
    } else NULL
    out <- bind_rows(rows1, rows2)
  } else {
    cat(sprintf("  [skip %-22s schema not recognised]\n", sheet))
    return(NULL)
  }

  out |>
    filter(!is.na(transition_id), nzchar(trimws(transition_id)),
           !is.na(precursor_mz), !is.na(product_mz))
}

# ---- 4. read all panel sheets ----------------------------------------------
panel_libs <- excel_sheets(PANELS) |>
  set_names() |>
  map(function(s) {
    out <- tryCatch(normalize_panel(s),
                    error = function(e) { cat("  [error", s, "]:", e$message, "\n"); NULL })
    if (!is.null(out) && nrow(out) > 0)
      cat(sprintf("  %-22s %4d rows\n", s, nrow(out)))
    out
  }) |>
  compact()

panels_combined <- bind_rows(panel_libs)
cat(sprintf("[panels total]  %d rows\n\n", nrow(panels_combined)))

# ---- 5. combine base + panels, dedup with conflict flag --------------------
all_lib <- bind_rows(base_lib, panels_combined) |>
  mutate(
    transition_id = trimws(as.character(transition_id)),
    Polarity      = trimws(as.character(Polarity)),
    precursor_mz  = round(suppressWarnings(as.numeric(precursor_mz)), 4),
    product_mz    = round(suppressWarnings(as.numeric(product_mz)), 4)
  )

cat(sprintf("[combined raw]  %d rows  (before dedup)\n", nrow(all_lib)))

dedup_key <- c("transition_id", "precursor_mz", "product_mz", "Polarity")

dedup <- all_lib |>
  group_by(across(all_of(dedup_key))) |>
  summarise(
    # IMPORTANT: compute conflict counts BEFORE re-using the column names
    # below. dplyr summarise evaluates top-to-bottom and references the
    # most recent binding of a name.
    n_distinct_CE   = n_distinct(collision_eV[!is.na(collision_eV)]),
    n_distinct_cone = n_distinct(cone_V[!is.na(cone_V)]),
    n_distinct_type = n_distinct(Type[!is.na(Type)]),
    n_sources       = n(),
    sources         = paste(unique(source), collapse = " | "),
    collision_eV    = first(stats::na.omit(collision_eV)) %||% NA_real_,
    cone_V          = first(stats::na.omit(cone_V))       %||% NA_real_,
    Type            = first(stats::na.omit(Type))         %||% NA_character_,
    Lipid           = first(stats::na.omit(Lipid))        %||% NA_character_,
    .groups = "drop"
  ) |>
  mutate(
    has_conflict = n_distinct_CE > 1 | n_distinct_cone > 1 | n_distinct_type > 1
  ) |>
  select(-n_sources, -n_distinct_CE, -n_distinct_cone, -n_distinct_type)

cat(sprintf("[dedup]         %d rows  (will flag after Lipid fallback)\n\n",
            nrow(dedup)))

# ---- 6. standardize Lipid → Class promoted, Lipid coarsened ----------------
lipid_map <- c(
  # Polar
  "Polar_neg" = "Polar", "Polar_pos" = "Polar", "Polar_met" = "Polar",
  "Polar metabolite" = "Polar", "Polar Pos" = "Polar", "Polar" = "Polar",
  # Oxylipin
  "Oxylipin" = "Oxylipin", "oxylipin" = "Oxylipin",
  "SPM" = "Oxylipin", "SPMs" = "Oxylipin", "cysLT" = "Oxylipin",
  # Sphingolipid
  "SL" = "Sphingolipid", "SL_pos" = "Sphingolipid",
  "Sphingolipid" = "Sphingolipid", "Sphingolipids" = "Sphingolipid",
  "Ceramides" = "Sphingolipid",
  # Phospholipid
  "PL" = "Phospholipid", "PL_pos" = "Phospholipid",
  "Phospholipid" = "Phospholipid", "Phospholipids" = "Phospholipid",
  "PC" = "Phospholipid", "PE" = "Phospholipid", "PI" = "Phospholipid",
  "PG" = "Phospholipid", "PS" = "Phospholipid",
  # Lysophospholipid
  "LP" = "Lysophospholipid", "LPC" = "Lysophospholipid",
  "LPI" = "Lysophospholipid", "LPG" = "Lysophospholipid",
  "LPS" = "Lysophospholipid", "LPE" = "Lysophospholipid",
  "Lysophospholipids" = "Lysophospholipid",
  # FFA / PUFA
  "FFA" = "FFA", "PUFA" = "FFA",
  # MAG/DAG/TAG
  "MAG/DAG/TAG" = "MAG_DAG_TAG", "MAG" = "MAG_DAG_TAG",
  # Cannabinoid
  "Cannabinoids" = "Cannabinoid", "Cannabinoid" = "Cannabinoid",
  # Cholesterol ester
  "Cholesterol ester"  = "Cholesterol_ester",
  "Cholesterol esters" = "Cholesterol_ester",
  "Cholesterol_ester"  = "Cholesterol_ester",
  # Oxysterol
  "Oxysterol" = "Oxysterol",
  # Xenobiotic / IS / Technical / Ink / Marco
  "Xenobiotic" = "Xenobiotic", "PFAS" = "Xenobiotic",
  "IS" = "IS", "Ink" = "Technical", "Technical" = "Technical",
  "Marco_panel" = "Other", "Neuro" = "Other",
  "Lipidomics" = "Lipid", "Lipid" = "Lipid"
)

# Sheet-name → Lipid fallback for panels that lack a Lipid_class column.
# Applied only when Lipid is NA/blank.
sheet_to_lipid <- c(
  "technical"           = "Technical",
  "Oxylipins"           = "Oxylipin",
  "resolvins"           = "Oxylipin",
  "cysLTs"              = "Oxylipin",
  "oxylipins_TqA_optim" = "Oxylipin",
  "Fatty Acids"         = "FFA",
  "Cholesterol esters"  = "Cholesterol_ester",
  "Lysophospholipids"   = "Lysophospholipid",
  "PA_LPA"              = "Phospholipid",
  "Phospholipids"       = "Phospholipid",
  "DAG_TAG"             = "MAG_DAG_TAG",
  "Cannabinoids"        = "Cannabinoid",
  "PolarPos"            = "Polar",
  "Xenobiotics"         = "Xenobiotic",
  "PolarBasicNeg"       = "Polar",
  "Sphingo_ceramides"   = "Sphingolipid"
)

infer_from_sources <- function(srcs) {
  # `srcs` is a "panelA | panelB | ..." string; pick the first known mapping
  toks <- trimws(unlist(strsplit(srcs, "\\s*\\|\\s*")))
  hits <- sheet_to_lipid[toks]
  hits <- hits[!is.na(hits)]
  if (length(hits) == 0) NA_character_ else unname(hits[1])
}

# Capture original specific value as Class, then coarsen Lipid via the map.
# When Lipid is NA/blank, fall back to the source sheet name.
dedup <- dedup |>
  mutate(
    Class = Lipid,
    Lipid = ifelse(!is.na(Lipid) & nzchar(Lipid) & Lipid %in% names(lipid_map),
                   unname(lipid_map[Lipid]),
                   Lipid),
    Lipid = ifelse(is.na(Lipid) | !nzchar(Lipid),
                   vapply(sources, infer_from_sources, character(1)),
                   Lipid)
  )

# ---- 7. read metadata sheets ----------------------------------------------
read_meta_safe <- function(sheet, cols) {
  d <- suppressWarnings(suppressMessages(read_excel(META, sheet)))
  missing <- setdiff(cols, names(d))
  for (m in missing) d[[m]] <- NA
  d |> select(all_of(cols))
}

# NEW oxylipin_metadata
oxy_cols <- c("Processing_name", "PUFA", "Primary",
              "Enzymatic_pathway", "Immediate_Pathway",
              "Autoxidation", "Enzyme_1", "Precursors", "Omega_6", "Omega_3")
oxy_meta <- read_meta_safe("NEW oxylipin_metadata", oxy_cols) |>
  rename(meta_name = Processing_name) |>
  filter(!is.na(meta_name) & nzchar(meta_name))

# SL_metadata
sl_cols <- c("Processing_name", "Class_abbrev", "DeNovo_synthesis",
             "Enzyme_list", "GBA-GALC", "PPAP2A", "ENPP7-SMPD",
             "DEGS", "CERK", "CERS1", "CERS2", "CERS3", "CERS4",
             "CERS5", "CERS6", "UGT8-UGCG", "GLA-NEU", "GLB-NEU",
             "B4GALT", "SGMS", "SGMS1", "SPHK", "KDSR", "SPTLC",
             "SGPP", "ACER", "Unkown")
sl_meta <- read_meta_safe("SL_metadata", sl_cols) |>
  rename(meta_name = Processing_name) |>
  filter(!is.na(meta_name) & nzchar(meta_name))

# Polar metadata (combine Neg + Pos)
polar_cols <- c("Processing_name", "Super_Pathway", "Sub_Pathway")
polar_meta <- bind_rows(
    read_meta_safe("PolarNeg_metadata", polar_cols),
    read_meta_safe("PolarPos_metadata", polar_cols)
  ) |>
  rename(meta_name = Processing_name) |>
  filter(!is.na(meta_name) & nzchar(meta_name)) |>
  distinct(meta_name, .keep_all = TRUE)

cat(sprintf("[metadata read but DROPPED from output]  oxy=%d  sl=%d  polar=%d\n",
            nrow(oxy_meta), nrow(sl_meta), nrow(polar_meta)))

# ---- 7b. Class from LC-MSMS feature_metadata sheets -----------------------
read_class_from_lcms <- function(f) {
  if (!file.exists(f)) {
    cat(sprintf("  [skip lcms %s -- not found]\n", basename(f)))
    return(NULL)
  }
  d <- suppressWarnings(suppressMessages(read_excel(f, "feature_metadata")))
  if (!"class" %in% names(d) || !"Processing_name" %in% names(d)) return(NULL)
  tibble(meta_name  = as.character(d$Processing_name),
         lcms_class = as.character(d$class)) |>
    filter(!is.na(meta_name) & nzchar(meta_name) &
           !is.na(lcms_class) & nzchar(lcms_class))
}

lcms_class <- bind_rows(lapply(LCMS_FILES, read_class_from_lcms)) |>
  distinct(meta_name, .keep_all = TRUE)
cat(sprintf("[lcms class]  %d entries from %d LC-MSMS files\n",
            nrow(lcms_class), length(LCMS_FILES)))

# ---- 8. join metadata by name (split on ' || ') ----------------------------
# Each transition_id may be "A || B || C"; we explode to long, join, and pick
# the first non-NA hit per transition.
long_keys <- dedup |>
  mutate(.row = row_number(),
         .toks = str_split(transition_id, "\\s*\\|\\|\\s*")) |>
  select(.row, .toks) |>
  unnest(.toks) |>
  mutate(.toks = trimws(.toks)) |>
  filter(nzchar(.toks))

attach_meta <- function(meta_df) {
  long_keys |>
    left_join(meta_df, by = c(".toks" = "meta_name")) |>
    group_by(.row) |>
    summarise(across(-c(.toks),
                     ~ { x <- .x[!is.na(.x) & nzchar(as.character(.x))]
                         if (length(x) > 0) x[1] else NA }),
              .groups = "drop")
}

# Only LC-MSMS class is joined (used to refine `Class`).
# All other metadata sheets are read for reference but NOT merged into output.
lcms_join  <- attach_meta(lcms_class)

dedup <- dedup |>
  mutate(.row = row_number()) |>
  left_join(lcms_join,  by = ".row") |>
  # Promote LC-MSMS class over panel-derived Class when available
  mutate(Class = ifelse(!is.na(lcms_class) & nzchar(lcms_class),
                        lcms_class, Class)) |>
  select(-.row, -lcms_class)

# ---- 8b. compute final flag (after Lipid fallback + LC-MSMS Class merge) --
dedup <- dedup |>
  mutate(
    has_blank = is.na(collision_eV) | is.na(cone_V) |
                is.na(Polarity) | !nzchar(as.character(Polarity)) |
                is.na(Type)     | !nzchar(as.character(Type)) |
                is.na(Lipid)    | !nzchar(as.character(Lipid)),
    flag = case_when(
      has_conflict & has_blank ~ "conflict_blank",
      has_conflict             ~ "conflict",
      has_blank                ~ "blank",
      TRUE                     ~ "ok"
    )
  ) |>
  select(-has_blank, -has_conflict)

cat(sprintf("[flag]  ok=%d  conflict=%d  blank=%d  conflict_blank=%d\n",
            sum(dedup$flag == "ok"),
            sum(dedup$flag == "conflict"),
            sum(dedup$flag == "blank"),
            sum(dedup$flag == "conflict_blank")))

# ---- 9. final column order -------------------------------------------------
# Lean schema: core MS parameters + LC-MSMS-refined Class + provenance.
# Detailed metadata (oxylipin/SL/polar) intentionally omitted — use the
# search_ion_library.R tool to look transitions up.
final <- dedup |>
  select(any_of(c(
    "transition_id", "precursor_mz", "product_mz",
    "collision_eV", "cone_V", "Polarity",
    "Type", "Lipid", "Class",
    "sources", "flag"
  )))

# ---- 10. write -------------------------------------------------------------
write_csv(final, OUT_CSV, na = "")
cat(sprintf("\n[written]  %s\n           %d rows  x  %d cols\n",
            OUT_CSV, nrow(final), ncol(final)))

# Per-flag summary at the end so it's easy to scan
cat("\n[flag summary]\n")
print(final |> count(flag))

cat("\n[Lipid (coarse) summary]\n")
print(final |> count(Lipid, sort = TRUE))
