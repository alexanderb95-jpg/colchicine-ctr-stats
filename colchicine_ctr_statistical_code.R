# STUDY-21-001175 — full statistical analysis code (extracted from CTR R Markdown)
# Source: CTR_Submission_Results_Protocol2101175.Rmd
# Data: ProtocolNo2101175PRM_DATA_2026-05-18_2135.csv
# Helpers: protocol2101175_redcap_helpers.R

# =============================================================================
# protocol2101175_redcap_helpers.R
# =============================================================================
# Helpers for Protocol 2101175 REDCap CTR report (labels, therapy class, tumor type).

#' Load CTCAE v5 term labels from REDCap companion .r export
load_ctcae_term_labels <- function(path = here::here("ProtocolNo2101175PRM_R_2026-05-18_2135.r")) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  if (length(lines) > 0 && startsWith(lines[[1]], "\ufeff")) {
    lines[[1]] <- sub("^\ufeff", "", lines[[1]])
  }
  start_ln <- which(grepl("^mapping_ae_term_ctcaev5\\s*=\\s*c\\(", lines))[1]
  end_ln <- which(grepl("^mapping_meddra_soc\\s*=\\s*c\\(", lines))[1] - 1
  if (is.na(start_ln) || is.na(end_ln)) {
    stop("Could not locate mapping_ae_term_ctcaev5 in ", path, call. = FALSE)
  }

  # Line-parse (companion file has a few non-standard keys that break source()).
  terms <- character()
  for (ln in lines[start_ln:end_ln]) {
    m <- regexec(
      '^[[:space:]]*"([^"]+)"[[:space:]]*=[[:space:]]*"([^"]*)"',
      ln,
      perl = TRUE
    )
    if (length(m) < 1L || m[[1]][[1L]] < 0L) next
    caps <- regmatches(ln, m)[[1]]
    if (length(caps) < 3 || is.na(caps[[3]])) next
    code <- caps[[2]]
    raw_val <- as.character(caps[[3]])
    if (!nzchar(raw_val)) next
    term <- trimws(strsplit(raw_val, "\t", fixed = TRUE)[[1]][1])
    if (!nzchar(term)) next
    terms[code] <- term
  }
  terms
}

#' Map one CTCAE code (numeric or character) to preferred term
ctcae_term_label <- function(code, term_map) {
  if (length(code) == 0) return(character(0))
  key <- as.character(code)
  out <- unname(term_map[key])
  miss <- is.na(out) & !is.na(code) & key != ""
  if (any(miss)) {
    out[miss] <- paste0("CTCAE code ", key[miss])
  }
  out
}

#' Extract screening laboratory values for a subject (first non-missing on screening visits).
protocol2101175_screening_labs <- function(dat, subject_id) {
  screening_events <- c("screening_28_days_arm_1", "screening_arm_2", "screening_arm_3")
  lab_cols <- c(
    "hemoglobin_g_dl", "creatinine_mg_dl", "ast_u_l", "alt_u_l",
    "bilirubin_total_mg_dl", "neutrophils", "c_reactive_protein_value"
  )

  screen_rows <- dat %>%
    dplyr::filter(
      .data$study_id == .env$subject_id,
      .data$redcap_event_name %in% screening_events
    )

  primary_row <- screen_rows %>%
    dplyr::filter(is.na(.data$redcap_repeat_instrument) | .data$redcap_repeat_instrument == "") %>%
    dplyr::slice(1)

  source_rows <- if (nrow(primary_row) > 0) primary_row else screen_rows %>% dplyr::slice(1)
  if (nrow(source_rows) == 0) {
    return(tibble::tibble())
  }

  source_rows %>%
    dplyr::select(dplyr::any_of(lab_cols)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "lab", values_to = "value") %>%
    dplyr::mutate(value = suppressWarnings(as.numeric(.data$value))) %>%
    dplyr::filter(!is.na(.data$value)) %>%
    dplyr::group_by(.data$lab) %>%
    dplyr::summarise(value = .data$value[[1]], .groups = "drop") %>%
    tidyr::pivot_wider(names_from = .data$lab, values_from = .data$value)
}

#' Identify organ-function values that fail protocol inclusion thresholds.
protocol2101175_organ_function_failures <- function(labs) {
  if (nrow(labs) == 0) return(character())

  fails <- character()
  if ("hemoglobin_g_dl" %in% names(labs) && !is.na(labs$hemoglobin_g_dl) && labs$hemoglobin_g_dl < 8) {
    fails <- c(fails, paste0("hemoglobin ", labs$hemoglobin_g_dl, " g/dL (required >=8 g/dL)"))
  }
  if ("neutrophils" %in% names(labs) && !is.na(labs$neutrophils) && labs$neutrophils < 1.5) {
    fails <- c(fails, paste0("ANC ", labs$neutrophils, " x10^3/uL (required >=1.5 x10^3/uL)"))
  }
  if ("creatinine_mg_dl" %in% names(labs) && !is.na(labs$creatinine_mg_dl) && labs$creatinine_mg_dl > 1.5) {
    fails <- c(fails, paste0("creatinine ", labs$creatinine_mg_dl, " mg/dL (required <=1.5 mg/dL)"))
  }
  fails
}

#' Build publication-ready screen-failure reason/detail from REDCap screening row.
protocol2101175_screen_fail_labels <- function(row, labs = NULL) {
  reason_labels <- c(
    reason_for_screen_fail___1 = "Inclusion criteria",
    reason_for_screen_fail___2 = "Exclusion criteria",
    reason_for_screen_fail___3 = "Surpassed screening window",
    reason_for_screen_fail___4 = "Physician decision to discontinue",
    reason_for_screen_fail___5 = "Patient decision to discontinue",
    reason_for_screen_fail___6 = "Other"
  )

  inc_fail_labels <- c(
    "4" = "CRP below inclusion threshold",
    "8" = "Inadequate organ function"
  )

  exc_labels <- c(
    "1" = "Long-term colchicine use",
    "2" = "Active infection requiring systemic therapy",
    "3" = "Pregnant or breastfeeding",
    "4" = "Prior cancer treatment within 30 days / 5 half-lives",
    "5" = "Prior/concurrent malignancy interfering with assessment",
    "6" = "Active CNS metastases",
    "7" = "Investigational drug within 30 days",
    "8" = "Prohibited concomitant medications",
    "9" = "RA/vasculitis/SLE requiring active treatment",
    "10" = "Recent MI or NYHA class III+ heart failure"
  )

  checked_nums <- function(prefix) {
    cols <- grep(paste0("^", prefix, "___[0-9]+$"), names(row), value = TRUE)
    hit <- cols[!is.na(row[cols]) & row[cols] == 1]
    gsub(paste0(prefix, "___"), "", hit)
  }

  redcap_reason <- {
    hit <- names(reason_labels)[vapply(
      names(reason_labels),
      function(nm) !is.na(row[[nm]]) && row[[nm]] == 1,
      logical(1)
    )]
    if (length(hit) == 0) NA_character_ else paste(unname(reason_labels[hit]), collapse = "; ")
  }

  inc_nums <- checked_nums("inclusion_criteria_causing")
  exc_nums <- checked_nums("exclusion_criteria_causing")
  other <- row$if_other_please_mention_de[[1]]

  detail_parts <- character()
  clinical_reason <- NA_character_

  if ("8" %in% inc_nums) {
    clinical_reason <- "Inadequate organ function"
    organ_fails <- if (!is.null(labs)) protocol2101175_organ_function_failures(labs) else character()
    if (length(organ_fails) > 0) {
      detail_parts <- c(
        detail_parts,
        paste0(
          "Did not meet inclusion organ-function criteria; ",
          paste(organ_fails, collapse = "; ")
        )
      )
    } else {
      detail_parts <- c(
        detail_parts,
        "Did not meet inclusion organ-function criteria (specific abnormal value not available in export)"
      )
    }
  }

  if ("4" %in% inc_nums) {
    clinical_reason <- "CRP below inclusion threshold"
    crp_val <- if (!is.null(labs) && "c_reactive_protein_value" %in% names(labs)) {
      labs$c_reactive_protein_value
    } else {
      NA_real_
    }
    detail_parts <- c(
      detail_parts,
      if (!is.na(crp_val)) {
        paste0("CRP ", crp_val, " mg/L (required >5 mg/L)")
      } else {
        "Did not meet inclusion CRP criterion (CRP >5 mg/L)"
      }
    )
  }

  other_inc <- setdiff(inc_nums, c("4", "8"))
  if (length(other_inc) > 0) {
    if (is.na(clinical_reason)) clinical_reason <- "Did not meet inclusion criteria"
    detail_parts <- c(
      detail_parts,
      paste0("Additional inclusion criterion not met (REDCap code ", paste(other_inc, collapse = ", "), ")")
    )
  }

  if (length(exc_nums) > 0) {
    if (is.na(clinical_reason)) {
      clinical_reason <- "Met exclusion criteria"
    }
    detail_parts <- c(
      detail_parts,
      paste0(
        "Met exclusion criterion: ",
        paste(unname(exc_labels[exc_nums]), collapse = "; ")
      )
    )
  }

  if (!is.na(other) && nzchar(other)) {
    if (is.na(clinical_reason)) clinical_reason <- "Other"
    detail_parts <- c(detail_parts, other)
  }

  if (is.na(clinical_reason)) {
    clinical_reason <- redcap_reason
  }

  if (length(detail_parts) == 0) {
    detail_parts <- if (!is.na(redcap_reason)) {
      paste0("REDCap screen-fail category: ", redcap_reason)
    } else {
      NA_character_
    }
  }

  list(
    screen_fail_reason_redcap = redcap_reason,
    screen_fail_reason = clinical_reason,
    screen_fail_detail = paste(detail_parts, collapse = "; ")
  )
}

#' Map REDCap screening visit to protocol enrollment era.
protocol2101175_protocol_era <- function(screening_event) {
  dplyr::case_when(
    screening_event == "screening_28_days_arm_1" ~ "original",
    screening_event == "screening_arm_2" ~ "amended",
    screening_event == "screening_arm_3" ~ "cohort2",
    TRUE ~ NA_character_
  )
}

#' Human-readable label for protocol enrollment era.
protocol2101175_protocol_era_label <- function(era) {
  dplyr::case_when(
    era == "original" ~ "Original schedule (indefinite colchicine + surveillance imaging)",
    era == "amended" ~ "Amended schedule (fixed 2-week phase 0 intervention)",
    era == "cohort2" ~ "Cohort 2 (closed; not enrolled)",
    TRUE ~ NA_character_
  )
}

#' Enrolled patients under the original schedule (Kaplan–Meier–evaluable).
protocol2101175_km_eligible_ids <- function(disposition) {
  disposition %>%
    dplyr::filter(.data$enrolled, .data$protocol_era == "original") %>%
    dplyr::pull(.data$study_id) %>%
    sort()
}

#' Build accrual disposition (consent, enrollment, screen-fail) from screening visits.
build_protocol2101175_screen_disposition <- function(dat) {
  screening_events <- c("screening_28_days_arm_1", "screening_arm_2", "screening_arm_3")

  consent <- dat %>%
    dplyr::filter(!is.na(.data$date_informed_consent_sign)) %>%
    dplyr::group_by(.data$study_id) %>%
    dplyr::summarise(
      consent_date = min(.data$date_informed_consent_sign, na.rm = TRUE),
      .groups = "drop"
    )

  disp_rows <- dat %>%
    dplyr::filter(
      .data$redcap_event_name %in% screening_events,
      !is.na(.data$was_the_subject_enrolled)
    ) %>%
    dplyr::group_by(.data$study_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  if (nrow(disp_rows) == 0) {
    stop("No screening disposition rows with was_the_subject_enrolled found.", call. = FALSE)
  }

  label_mat <- lapply(seq_len(nrow(disp_rows)), function(i) {
    row <- disp_rows[i, , drop = FALSE]
    sid <- row$study_id[[1]]
    labs <- protocol2101175_screening_labs(dat, sid)
    protocol2101175_screen_fail_labels(row, labs)
  })

  consent %>%
    dplyr::left_join(
      disp_rows %>%
        dplyr::transmute(
          study_id = .data$study_id,
          screening_redcap_event = .data$redcap_event_name,
          enrolled_flag = .data$was_the_subject_enrolled,
          enroll_date = .data$date_of_enrollment
        ),
      by = "study_id"
    ) %>%
    dplyr::mutate(
      protocol_era = protocol2101175_protocol_era(.data$screening_redcap_event),
      protocol_era_label = protocol2101175_protocol_era_label(.data$protocol_era),
      screen_fail_reason_redcap = vapply(label_mat, `[[`, character(1), "screen_fail_reason_redcap"),
      screen_fail_reason = vapply(label_mat, `[[`, character(1), "screen_fail_reason"),
      screen_fail_detail = vapply(label_mat, `[[`, character(1), "screen_fail_detail"),
      enrolled = .data$enrolled_flag == 1,
      screen_fail = .data$enrolled_flag == 0,
      km_evaluable = .data$enrolled_flag == 1 & .data$protocol_era == "original",
      disposition = dplyr::case_when(
        .data$enrolled ~ "Enrolled",
        .data$screen_fail ~ "Screen failure",
        TRUE ~ "Unknown"
      )
    ) %>%
    dplyr::arrange(.data$study_id)
}

#' Summarize screen-failure reasons for accrual table
summarise_screen_fail_reasons <- function(disposition) {
  sf <- disposition %>% dplyr::filter(.data$screen_fail)
  if (nrow(sf) == 0) {
    return(tibble::tibble(
      `Screen-failure reason` = character(),
      `N` = integer(),
      `Details` = character()
    ))
  }

  sf %>%
    dplyr::group_by(.data$screen_fail_reason) %>%
    dplyr::summarise(
      `N` = dplyr::n(),
      `Details` = paste(unique(stats::na.omit(.data$screen_fail_detail)), collapse = " | "),
      .groups = "drop"
    ) %>%
    dplyr::rename(`Screen-failure reason` = .data$screen_fail_reason) %>%
    dplyr::arrange(dplyr::desc(.data$`N`))
}

#' Classify a prior-therapy row into a publication-friendly class
classify_prior_therapy <- function(therapy_name) {
  t <- tolower(trimws(therapy_name))
  dplyr::case_when(
    is.na(t) | t == "" ~ NA_character_,
    stringr::str_detect(t, "radiotherapy|^rt |^xrt|sbrt|rx dose|neoadjuvant rt|adjuvant rt") ~ "Radiation therapy",
    stringr::str_detect(
      t,
      "resection|ectomy|nephroureterectomy|exenteration|cystectomy|omentectomy|laparoscopic|excision"
    ) ~ "Surgery",
    stringr::str_detect(t, "pembrolizumab|nivolumab|atezolizumab|ipilimumab|durvalumab") ~ "Immunotherapy",
    stringr::str_detect(t, "bevacizumab|panitumumab|enfortumab|rucaparib|niraparib|olaparib|doxil") ~ "Targeted therapy",
    stringr::str_detect(
      t,
      "carboplatin|cisplatin|paclitaxel|taxol|taxotere|gemcitabine|oxaliplatin|irinotecan|fluorouracil|xeloda|capecitabine|doxorubicin|etoposide|topotecan|cytoxan|folinic"
    ) ~ "Chemotherapy",
    TRUE ~ "Other"
  )
}

#' RECIST target-lesion codes → short labels (REDCap mapping_target_lesions)
recist_response_label <- function(code) {
  mapping <- c(
    "1" = "CR",
    "2" = "PR",
    "3" = "SD",
    "4" = "PD",
    "5" = "Non-evaluable"
  )
  out <- unname(mapping[as.character(code)])
  out[is.na(out) & !is.na(code) & as.character(code) != ""] <- paste0("Code ", code)
  out
}

#' Build per-patient PFS/OS from REDCap export (RECIST BOR, discontinuation, final status).
#'
#' The dedicated progression_free_survival_pfs_overall_survival_os form was not completed
#' in this export; endpoints are derived from available tumor assessments and status fields.
build_protocol2101175_survival <- function(dat, enrolled_ids) {
  d <- dat %>% dplyr::filter(study_id %in% enrolled_ids)

  index <- d %>%
    dplyr::filter(!is.na(.data$date_of_enrollment)) %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(index_date = min(.data$date_of_enrollment, na.rm = TRUE), .groups = "drop")

  ip_start <- d %>%
    dplyr::filter(!is.na(.data$date_ip_dispensed)) %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(ip_start = min(.data$date_ip_dispensed, na.rm = TRUE), .groups = "drop")

  pd_bor <- d %>%
    dplyr::filter(
      .data$redcap_repeat_instrument == "best_overall_response",
      .data$target_lesions == 4
    ) %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(pd_bor = min(.data$date_of_overall_response, na.rm = TRUE), .groups = "drop")

  pd_disc <- d %>%
    dplyr::filter(
      .data$reason_for_discontinuation_from_protocol_therapy___2 == 1 |
        .data$reason_for_removal_of_subj___4 == 1
    ) %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(
      pd_disc = min(
        dplyr::coalesce(
          .data$discontinuation_from_protocol_therapy,
          .data$date_of_discontinuation_from_protocol_activities
        ),
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  death <- d %>%
    dplyr::filter(.data$subject_status_2 == 1, !is.na(.data$date_of_final_status)) %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(death_date = min(.data$date_of_final_status, na.rm = TRUE), .groups = "drop")

  last_alive <- d %>%
    dplyr::filter(!is.na(.data$date_subject_last_contacte)) %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(
      last_alive = max(.data$date_subject_last_contacte, na.rm = TRUE),
      .groups = "drop"
    )

  fu_date_cols <- c(
    "date_of_overall_response", "date_of_target_lesions_eva", "date_c_reactive_protein_co",
    "date_physical_exam_perform", "date_of_discontinuation_from_protocol_activities",
    "discontinuation_from_protocol_therapy", "end_date_monitoring",
    "date_subject_last_contacte", "date_of_final_status"
  )

  fu <- d %>%
    dplyr::select(study_id, dplyr::any_of(fu_date_cols)) %>%
    tidyr::pivot_longer(-study_id, values_to = "dt") %>%
    dplyr::mutate(dt = suppressWarnings(lubridate::ymd(.data$dt))) %>%
    dplyr::filter(!is.na(.data$dt)) %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(last_fu = max(.data$dt, na.rm = TRUE), .groups = "drop")

  bor <- d %>%
    dplyr::filter(
      .data$redcap_repeat_instrument == "best_overall_response",
      !is.na(.data$target_lesions)
    ) %>%
    dplyr::mutate(
      best_resp = recist_response_label(.data$target_lesions),
      bor_dt = dplyr::coalesce(.data$date_of_overall_response, .data$date_of_target_lesions_eva)
    ) %>%
    dplyr::arrange(.data$study_id, .data$bor_dt) %>%
    dplyr::group_by(study_id) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(study_id, best_recist = best_resp, last_bor_date = bor_dt)

  status <- d %>%
    dplyr::filter(.data$subject_status_2 != "") %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(
      final_status_code = .data$subject_status_2[which.max(dplyr::coalesce(
        .data$date_of_final_status,
        as.Date("1970-01-01")
      ))],
      .groups = "drop"
    )

  status_lab <- c("0" = "Alive", "1" = "Dead", "9" = "Lost to follow-up")

  index %>%
    dplyr::left_join(ip_start, by = "study_id") %>%
    dplyr::left_join(pd_bor, by = "study_id") %>%
    dplyr::left_join(pd_disc, by = "study_id") %>%
    dplyr::left_join(death, by = "study_id") %>%
    dplyr::left_join(last_alive, by = "study_id") %>%
    dplyr::left_join(fu, by = "study_id") %>%
    dplyr::left_join(bor, by = "study_id") %>%
    dplyr::left_join(status, by = "study_id") %>%
    dplyr::mutate(
      time_zero = dplyr::coalesce(.data$ip_start, .data$index_date),
      progression_date = pmin(.data$pd_bor, .data$pd_disc, na.rm = TRUE),
      pfs_event = !is.na(.data$progression_date),
      pfs_end_date = dplyr::if_else(
        .data$pfs_event,
        .data$progression_date,
        dplyr::coalesce(.data$last_bor_date, .data$last_fu, .data$last_alive)
      ),
      pfs_days = as.numeric(.data$pfs_end_date - .data$time_zero),
      pfs_months = .data$pfs_days / 30.44,
      os_event = !is.na(.data$death_date),
      os_end_date = dplyr::if_else(
        .data$os_event,
        .data$death_date,
        dplyr::coalesce(.data$last_alive, .data$last_fu)
      ),
      os_days = as.numeric(.data$os_end_date - .data$time_zero),
      os_months = .data$os_days / 30.44,
      final_status = unname(status_lab[as.character(.data$final_status_code)])
    ) %>%
    dplyr::filter(!is.na(.data$time_zero), !is.na(.data$pfs_end_date), .data$pfs_days >= 0) %>%
    dplyr::arrange(.data$study_id)
}

#' Best % change in target-lesion sum of diameters (SOD) for enrolled patients.
#'
#' Scan sets are inferred within target_tumor_lesion_details by resets of lesion_number_t
#' to 1. First set = baseline; subsequent sets = post-baseline. Best % change is the
#' minimum post-baseline change (most negative = maximum shrinkage). Patients with
#' fewer than 2 scan sets are excluded.
build_protocol2101175_tumor_change <- function(dat, enrolled_ids) {
  les <- dat %>%
    dplyr::filter(
      .data$study_id %in% enrolled_ids,
      .data$redcap_repeat_instrument == "target_tumor_lesion_details"
    ) %>%
    dplyr::mutate(
      lesion_number_t = suppressWarnings(as.numeric(.data$lesion_number_t)),
      ld = suppressWarnings(as.numeric(.data$longest_diameter_t))
    ) %>%
    dplyr::filter(!is.na(.data$ld), !is.na(.data$lesion_number_t)) %>%
    dplyr::arrange(.data$study_id, .data$redcap_repeat_instance)

  if (nrow(les) == 0) {
    return(tibble::tibble(
      study_id = character(),
      baseline_sod = numeric(),
      best_pct_change = numeric(),
      max_shrinkage_pct = numeric(),
      n_post_scans = integer()
    ))
  }

  assign_scan_sets <- function(df) {
    ln <- df$lesion_number_t
    set_id <- integer(nrow(df))
    cur <- 0L
    prev <- NA_real_
    for (i in seq_len(nrow(df))) {
      if (is.na(prev)) {
        cur <- 1L
      } else if (ln[[i]] == 1 && prev != 1) {
        cur <- cur + 1L
      }
      set_id[[i]] <- cur
      prev <- ln[[i]]
    }
    df$scan_set <- set_id
    df
  }

  les %>%
    split(.$study_id) %>%
    lapply(assign_scan_sets) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(study_id = as.character(.data$study_id)) %>%
    dplyr::group_by(.data$study_id, .data$scan_set) %>%
    dplyr::summarise(sod = sum(.data$ld), .groups = "drop") %>%
    dplyr::group_by(.data$study_id) %>%
    dplyr::filter(dplyr::n_distinct(.data$scan_set) >= 2) %>%
    dplyr::arrange(.data$scan_set, .by_group = TRUE) %>%
    dplyr::mutate(
      baseline_sod = .data$sod[.data$scan_set == min(.data$scan_set)][[1]],
      pct_change = dplyr::if_else(
        .data$scan_set == min(.data$scan_set),
        NA_real_,
        100 * (.data$sod - .data$baseline_sod) / .data$baseline_sod
      )
    ) %>%
    dplyr::filter(!is.na(.data$pct_change)) %>%
    dplyr::summarise(
      baseline_sod = dplyr::first(.data$baseline_sod),
      best_pct_change = min(.data$pct_change, na.rm = TRUE),
      max_shrinkage_pct = pmax(-min(.data$pct_change, na.rm = TRUE), 0),
      n_post_scans = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(study_id = as.character(.data$study_id)) %>%
    dplyr::arrange(.data$study_id)
}

#' Infer primary tumor type from prior therapy / surgery text (no dedicated REDCap diagnosis field)
infer_primary_tumor_type <- function(study_id, therapy_names) {
  t <- tolower(paste(unique(therapy_names[!is.na(therapy_names) & therapy_names != ""]), collapse = " | "))
  if (t == "") return(NA_character_)

  out <- dplyr::case_when(
    stringr::str_detect(t, "ovarian") ~ "Ovarian cancer",
    stringr::str_detect(t, "olaparib|rucaparib|niraparib") ~ "Ovarian cancer",
    stringr::str_detect(t, "enfortumab|nephroureterectomy|cystoprostatectomy") ~ "Urothelial carcinoma",
    stringr::str_detect(t, "rectum|abdominoperineal") ~ "Rectal adenocarcinoma",
    stringr::str_detect(t, "oral cavity|skull base") ~ "Head and neck squamous cell carcinoma",
    stringr::str_detect(t, "vulva|pelvic exenteration") ~ "Cervical cancer",
    stringr::str_detect(t, "anal") ~ "Anal squamous cell carcinoma",
    stringr::str_detect(t, "topotecan|gemcitabine") &
      stringr::str_detect(t, "nivolumab|pembrolizumab|ipilimumab") ~ "Urothelial carcinoma",
    stringr::str_detect(t, "carboplatin") &
      stringr::str_detect(t, "taxol|paclitaxel") &
      stringr::str_detect(t, "bevacizumab|doxil|atezolizumab") ~ "Ovarian cancer",
    stringr::str_detect(t, "etoposide") & stringr::str_detect(t, "carboplatin") ~ "Small cell carcinoma",
    TRUE ~ NA_character_
  )
  if (is.na(out)) {
    warning("Could not infer tumor type for study_id ", study_id, call. = FALSE)
    out <- "Not specified in REDCap export"
  }
  out
}

# =============================================================================
# R chunk: setup  (, include=FALSE)
# =============================================================================
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5)
set.seed(123)

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  library(knitr)
  library(kableExtra)
  library(survival)
  library(survminer)
})

# ---- Load REDCap export and label helpers ----
source(here::here("scripts", "protocol2101175_redcap_helpers.R"))
data_path <- here::here(
  "data", "clinical", "protocol_2101175",
  "ProtocolNo2101175PRM_DATA_2026-05-18_2135.csv"
)
raw <- read_csv(data_path, show_col_types = FALSE)
ctcae_term_map <- load_ctcae_term_labels()

# ---- Factor mappings (from REDCap companion script) ----
# Safe recode: never assign a full vector into is.na(out) slots (avoids label corruption).
map_chr <- function(x, mapping) {
  x_chr <- as.character(x)
  out <- unname(mapping[x_chr])
  unmatched <- is.na(out) & !is.na(x) & x_chr != ""
  if (any(unmatched)) {
    out[unmatched] <- paste0("Code ", x_chr[unmatched])
  }
  out[is.na(x) | x_chr == ""] <- NA_character_
  out
}

# Stop knit if categorical table cells do not sum to the stated total N.
assert_partition_sums <- function(table_name, component_counts, total_n) {
  comp <- as.integer(component_counts)
  if (any(is.na(comp))) {
    stop(table_name, ": missing component counts (NA).", call. = FALSE)
  }
  s <- sum(comp)
  if (s != total_n) {
    stop(
      sprintf(
        "%s: components (%s) sum to %d but total N = %d",
        table_name,
        paste(comp, collapse = " + "),
        s,
        total_n
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

cohort_lab <- c(
  "1" = "Cohort 1 low dose (0.6 mg BID)",
  "2" = "Cohort 1 high dose (0.6 mg TID)",
  "3" = "Cohort 2 (0.6 mg BID + nivolumab)"
)
gender_lab <- c("1" = "Male", "2" = "Female", "3" = "Other")
race_lab <- c(
  "1" = "Black", "2" = "White", "3" = "Asian", "4" = "American Indian",
  "5" = "Hawaiian/Pacific Islander", "6" = "Other/Unknown"
)
ethnicity_lab <- c(
  "1" = "Hispanic or Latino", "2" = "NOT Hispanic or Latino", "3" = "Unknown / Not Reported"
)
ecog_lab <- c(
  "0" = "0", "1" = "1", "2" = "2", "3" = "3", "4" = "4", "-99" = "Not collected"
)
enrolled_lab <- c("0" = "No", "1" = "Yes")
status_lab <- c("0" = "Alive", "1" = "Dead", "9" = "Lost to follow-up")
ae_grade_lab <- c(
  "1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5"
)
sae_lab <- c("0" = "No", "1" = "Yes")
attr_lab <- c(
  "1" = "Unrelated", "2" = "Unlikely", "3" = "Possible", "4" = "Probable", "5" = "Definite", "6" = "N/A"
)

dat <- raw %>%
  mutate(
    study_id = as.integer(study_id),
    cohort_label = map_chr(cohort, cohort_lab),
    gender_label = map_chr(gender, gender_lab),
    race_label = map_chr(race, race_lab),
    ethnicity_label = map_chr(ethnicity, ethnicity_lab),
    enrolled_label = map_chr(was_the_subject_enrolled, enrolled_lab),
    ecog_label = map_chr(ecog_performance_status, ecog_lab),
    subject_status_label = map_chr(subject_status_2, status_lab),
    ae_grade_label = map_chr(aegrade, ae_grade_lab),
    sae_label = map_chr(is_this_a_sae, sae_lab),
    ae_attribution_label = map_chr(ae_attribution_colchicine, attr_lab),
    crp_value = suppressWarnings(as.numeric(c_reactive_protein_value)),
    ferritin_value_num = suppressWarnings(as.numeric(ferritin_value)),
    esr_value_num = suppressWarnings(as.numeric(erythrocyte_sedimentation_value)),
    date_of_birth_dt = suppressWarnings(ymd(date_of_birth)),
    date_enrollment_dt = suppressWarnings(ymd(date_of_enrollment)),
    date_final_status_dt = suppressWarnings(ymd(date_of_final_status)),
    treatment_start_dt = suppressWarnings(ymd(treatment_start_date))
  )

# ---- Cohort / visit helpers ----
baseline_events <- c(
  "cycle_1_day_1__3_d_arm_1", "cycle_1_day_1_arm_2", "cycle_1_day_1_arm_3"
)
post_baseline_events <- c(
  "cycle_1_day_15__3_arm_1",
  "cycle_1_day_8_arm_2", "cycle_1_day_15_arm_2",
  "cycle_1_day_8_arm_3", "cycle_1_day_15_arm_3"
)
screening_events <- c("screening_28_days_arm_1", "screening_arm_2", "screening_arm_3")

screen_disposition <- build_protocol2101175_screen_disposition(dat)
consented_ids <- screen_disposition %>%
  filter(!is.na(consent_date)) %>%
  pull(study_id) %>%
  sort()
enrolled_ids <- screen_disposition %>%
  filter(enrolled) %>%
  pull(study_id) %>%
  sort()
screen_fail_ids <- screen_disposition %>%
  filter(screen_fail) %>%
  pull(study_id) %>%
  sort()

enrolled_original_ids <- screen_disposition %>%
  filter(enrolled, protocol_era == "original") %>%
  pull(study_id) %>%
  sort()
enrolled_amended_ids <- screen_disposition %>%
  filter(enrolled, protocol_era == "amended") %>%
  pull(study_id) %>%
  sort()
km_eligible_ids <- protocol2101175_km_eligible_ids(screen_disposition)

consented_n <- length(consented_ids)
enrolled_n <- length(enrolled_ids)
enrolled_original_n <- length(enrolled_original_ids)
enrolled_amended_n <- length(enrolled_amended_ids)
km_eligible_n <- length(km_eligible_ids)
screen_fail_n <- length(screen_fail_ids)
screen_fail_original_n <- screen_disposition %>%
  filter(screen_fail, protocol_era == "original") %>%
  nrow()
screen_fail_amended_n <- screen_disposition %>%
  filter(screen_fail, protocol_era == "amended") %>%
  nrow()
screened_n <- consented_n

assert_partition_sums(
  "Accrual disposition (consented vs enrolled + screen-fail)",
  c(enrolled_n, screen_fail_n),
  consented_n
)
assert_partition_sums(
  "Enrolled patients by protocol schedule",
  c(enrolled_original_n, enrolled_amended_n),
  enrolled_n
)
assert_partition_sums(
  "Screen failures by protocol schedule",
  c(screen_fail_original_n, screen_fail_amended_n),
  screen_fail_n
)

screen_fail_summary <- summarise_screen_fail_reasons(screen_disposition)

# ---- Primary endpoint: max % CRP decline ----
crp_long <- dat %>%
  filter(study_id %in% enrolled_ids, !is.na(crp_value)) %>%
  select(study_id, redcap_event_name, crp_value, cohort_label) %>%
  distinct()

derive_crp_endpoint <- function(id, d_long) {
  sub <- d_long %>% filter(study_id == id)
  base_row <- sub %>%
    filter(redcap_event_name %in% baseline_events) %>%
    arrange(redcap_event_name) %>%
    slice(1)
  if (nrow(base_row) == 0) {
    base_row <- sub %>%
      filter(redcap_event_name %in% screening_events) %>%
      arrange(redcap_event_name) %>%
      slice(1)
  }
  if (nrow(base_row) == 0) return(NULL)
  baseline <- base_row$crp_value[1]
  if (is.na(baseline) || baseline <= 0) return(NULL)
  post <- sub %>%
    filter(redcap_event_name %in% post_baseline_events) %>%
    pull(crp_value)
  if (length(post) == 0) return(NULL)
  # Decline cannot be negative by definition; increases contribute 0% decline.
  pct_declines <- pmax(100 * (baseline - post) / baseline, 0)
  tibble(
    study_id = id,
    baseline_crp_mg_l = baseline,
    min_on_tx_crp_mg_l = min(post),
    max_pct_decline = max(pct_declines),
    n_crp_on_tx = length(post)
  )
}

crp_endpoint <- bind_rows(lapply(enrolled_ids, derive_crp_endpoint, crp_long))
evaluable_efficacy_n <- nrow(crp_endpoint)
evaluable_toxicity_n <- enrolled_n

# Primary one-sided t-test (protocol: alpha = 0.025 for Cohort 1 low dose)
primary_test <- if (evaluable_efficacy_n >= 2) {
  t.test(crp_endpoint$max_pct_decline, mu = 0, alternative = "greater")
} else {
  NULL
}

primary_mean <- mean(crp_endpoint$max_pct_decline, na.rm = TRUE)
primary_sd <- sd(crp_endpoint$max_pct_decline, na.rm = TRUE)
primary_p <- if (!is.null(primary_test)) primary_test$p.value else NA_real_
primary_t <- if (!is.null(primary_test)) primary_test$statistic else NA_real_
primary_ci_low <- if (!is.null(primary_test)) primary_test$conf.int[1] else NA_real_

cohens_d_primary <- if (evaluable_efficacy_n >= 2 && !is.na(primary_sd) && primary_sd > 0) {
  primary_mean / primary_sd
} else {
  NA_real_
}

posthoc_power_pct <- if (evaluable_efficacy_n >= 2 && !is.na(cohens_d_primary)) {
  n <- evaluable_efficacy_n
  ncp <- cohens_d_primary * sqrt(n)
  t_crit <- qt(1 - 0.025, df = n - 1)
  (1 - pt(t_crit, df = n - 1, ncp = ncp)) * 100
} else {
  NA_real_
}

# ---- Patient-level master (one row per enrolled subject) ----
# Demographics from screening visit; first non-missing value per field (not first table row).
patient_master <- tibble(study_id = enrolled_ids) %>%
  left_join(
    dat %>%
      filter(study_id %in% enrolled_ids, redcap_event_name %in% screening_events) %>%
      group_by(study_id) %>%
      summarise(
        gender_code = suppressWarnings(as.integer(first(na.omit(gender)))),
        race_code = suppressWarnings(as.integer(first(na.omit(race)))),
        ethnicity_code = suppressWarnings(as.integer(first(na.omit(ethnicity)))),
        ecog_code = suppressWarnings(as.integer(first(na.omit(ecog_performance_status)))),
        date_of_birth_dt = first(na.omit(date_of_birth_dt)),
        date_enrollment_dt = first(na.omit(date_enrollment_dt)),
        cohort_code = suppressWarnings(as.integer(first(na.omit(cohort)))),
        .groups = "drop"
      ),
    by = "study_id"
  ) %>%
  mutate(
    gender_label = map_chr(gender_code, gender_lab),
    race_label = map_chr(race_code, race_lab),
    ethnicity_label = map_chr(ethnicity_code, ethnicity_lab),
    ecog_label = map_chr(ecog_code, ecog_lab),
    cohort_label = map_chr(cohort_code, cohort_lab),
    age_years = if_else(
      !is.na(date_of_birth_dt) & !is.na(date_enrollment_dt),
      as.numeric(interval(date_of_birth_dt, date_enrollment_dt) / years(1)),
      NA_real_
    )
  )

if (nrow(patient_master) != enrolled_n) {
  stop("patient_master row count does not match enrolled_n.", call. = FALSE)
}

baseline <- patient_master

# Prior therapies → therapy class (repeat instrument)
prior_tx <- dat %>%
  filter(
    study_id %in% enrolled_ids,
    redcap_repeat_instrument == "prior_treatment_therapy",
    prior_treatment_therapies != ""
  ) %>%
  distinct(study_id, prior_treatment_therapies, therapy_start_date) %>%
  mutate(therapy_class = classify_prior_therapy(prior_treatment_therapies))

therapy_class_levels <- c(
  "Chemotherapy",
  "Immunotherapy",
  "Targeted therapy",
  "Radiation therapy",
  "Surgery",
  "Other"
)

therapy_class_levels_table1 <- c(
  "Chemotherapy",
  "Immunotherapy",
  "Targeted therapy",
  "Radiation therapy",
  "Surgery"
)

prior_tx_by_patient <- prior_tx %>%
  distinct(study_id, therapy_class) %>%
  filter(!is.na(therapy_class)) %>%
  mutate(received = TRUE) %>%
  pivot_wider(
    names_from = therapy_class,
    values_from = received,
    values_fill = list(received = FALSE)
  )

for (col in therapy_class_levels) {
  if (!col %in% names(prior_tx_by_patient)) {
    prior_tx_by_patient[[col]] <- FALSE
  }
}

prior_tx_summary <- prior_tx %>%
  group_by(therapy_class) %>%
  summarise(
    `N patients` = n_distinct(study_id),
    `N regimen entries` = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(`N patients`))

# Primary tumor type (inferred from prior therapy / surgery text in REDCap)
tumor_by_patient <- prior_tx %>%
  group_by(study_id) %>%
  summarise(
    primary_tumor_type = infer_primary_tumor_type(
      study_id[1],
      prior_treatment_therapies
    ),
    .groups = "drop"
  )

patient_master <- patient_master %>%
  left_join(tumor_by_patient, by = "study_id")

# Systemic therapy line count approximation:
# collapse same-date combinations into one regimen line and exclude RT/surgery/other.
prior_systemic_lines <- prior_tx %>%
  filter(therapy_class %in% c("Chemotherapy", "Immunotherapy", "Targeted therapy")) %>%
  mutate(systemic_start_date = suppressWarnings(ymd(therapy_start_date))) %>%
  group_by(study_id) %>%
  summarise(
    n_systemic_lines = if_else(
      any(!is.na(systemic_start_date)),
      n_distinct(systemic_start_date[!is.na(systemic_start_date)]),
      n_distinct(prior_treatment_therapies)
    ),
    .groups = "drop"
  ) %>%
  right_join(tibble(study_id = enrolled_ids), by = "study_id") %>%
  mutate(n_systemic_lines = replace_na(n_systemic_lines, 0L))

tumor_type_summary <- patient_master %>%
  count(primary_tumor_type, name = "N") %>%
  arrange(desc(N)) %>%
  mutate(
    `Cancer type (patient-level)` = paste0(primary_tumor_type, " (n = ", N, ")")
  )

assert_partition_sums(
  "Patient characteristics: primary cancer type",
  tumor_type_summary$N,
  enrolled_n
)

# ---- Adverse events ----
# TRAE-only: attribution to colchicine is Possible / Probable / Definite.
ae_all <- dat %>%
  filter(
    study_id %in% enrolled_ids,
    ae_term_ctcaev5 != "",
    ae_attribution_colchicine %in% c(3, 4, 5)
  ) %>%
  mutate(
    ae_grade_num = suppressWarnings(as.numeric(aegrade)),
    ae_term_display = ctcae_term_label(ae_term_ctcaev5, ctcae_term_map),
    ae_other_text = if ("if_you_entered_other_pleas" %in% names(.))
      trimws(as.character(if_you_entered_other_pleas))
    else
      NA_character_
  ) %>%
  mutate(
    ae_term_display = if_else(
      !is.na(ae_other_text) & ae_other_text != "" &
        (is.na(ae_term_display) |
          stringr::str_detect(tolower(ae_term_display), "^other|not elsewhere classified")),
      ae_other_text,
      ae_term_display
    )
  )

ae_by_term_grade <- ae_all %>%
  count(ae_term_display, ae_grade_label, name = "n") %>%
  arrange(desc(n))

# Worst CTCAE grade per patient per preferred term (do not count lower grades separately).
ae_worst_per_patient <- ae_all %>%
  filter(!is.na(ae_grade_num)) %>%
  group_by(study_id, ae_term_display) %>%
  slice_max(ae_grade_num, n = 1, with_ties = FALSE) %>%
  ungroup()

ae_by_term_grade_patient <- ae_worst_per_patient %>%
  count(ae_term_display, ae_grade_label, name = "N patients") %>%
  mutate(`N (%)` = paste0(`N patients`, " (", sprintf("%.1f", 100 * `N patients` / enrolled_n), "%)")) %>%
  arrange(desc(`N patients`), ae_term_display, suppressWarnings(as.numeric(ae_grade_label)))

ae_any_grade <- ae_all %>%
  distinct(study_id, ae_term_display) %>%
  count(study_id, name = "n_aes") %>%
  right_join(tibble(study_id = enrolled_ids), by = "study_id") %>%
  mutate(n_aes = replace_na(n_aes, 0L))

grade3plus_rate <- ae_worst_per_patient %>%
  filter(ae_grade_num >= 3) %>%
  distinct(study_id) %>%
  nrow() / enrolled_n * 100

trae_patient_n <- ae_all %>%
  distinct(study_id) %>%
  nrow()
trae_patient_pct <- 100 * trae_patient_n / enrolled_n

trae_diarrhea_n <- ae_worst_per_patient %>%
  filter(tolower(ae_term_display) == "diarrhea") %>%
  distinct(study_id) %>%
  nrow()
trae_diarrhea_pct <- 100 * trae_diarrhea_n / enrolled_n

# Helpers for inline text (also used in survival tables)
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  format(round(p, 3), nsmall = 3)
}
fmt_num <- function(x, digits = 1) format(round(x, digits), nsmall = digits)

.tbl_n <- 0L
.fig_n <- 0L
tbl_cap <- function(title) {
  .tbl_n <<- .tbl_n + 1L
  paste0("Table ", .tbl_n, ". ", title)
}
fig_cap <- function(title) {
  .fig_n <<- .fig_n + 1L
  paste0("Figure ", .fig_n, ". ", title)
}

# ---- Survival / status (PFS/OS derived from RECIST, discontinuation, final status) ----
# Restrict to original-schedule enrollments with serial imaging surveillance (Arm 1).
surv_dat <- build_protocol2101175_survival(dat, km_eligible_ids)

# Best % change in target-lesion SOD (exploratory; patients with ≥2 scan sets)
tumor_change <- build_protocol2101175_tumor_change(dat, enrolled_ids)

status_tbl <- surv_dat %>%
  transmute(
    study_id,
    status = coalesce(final_status, "Not recorded"),
    date_final = if_else(os_event, death_date, os_end_date)
  )

n_dead <- sum(surv_dat$os_event, na.rm = TRUE)
n_pfs_events <- sum(surv_dat$pfs_event, na.rm = TRUE)

surv_pfs <- surv_dat %>%
  filter(!is.na(pfs_months), pfs_months >= 0)

surv_os <- surv_dat %>%
  filter(!is.na(os_months), os_months >= 0)

fit_pfs <- if (nrow(surv_pfs) >= 3 && sum(surv_pfs$pfs_event) >= 1) {
  survfit(Surv(pfs_months, pfs_event) ~ 1, data = surv_pfs)
} else {
  NULL
}

fit_os <- if (nrow(surv_os) >= 3 && sum(surv_os$os_event) >= 1) {
  survfit(Surv(os_months, os_event) ~ 1, data = surv_os)
} else {
  NULL
}

median_surv_months <- function(fit) {
  if (is.null(fit)) return(list(median = NA_real_, lcl = NA_real_, ucl = NA_real_))
  s <- summary(fit)$table
  if (is.matrix(s)) {
    list(
      median = unname(s[1, "median"]),
      lcl = unname(s[1, "0.95LCL"]),
      ucl = unname(s[1, "0.95UCL"])
    )
  } else {
    list(median = NA_real_, lcl = NA_real_, ucl = NA_real_)
  }
}

pfs_med <- median_surv_months(fit_pfs)
os_med <- median_surv_months(fit_os)

pfs_os_tbl <- surv_dat %>%
  transmute(
    `Study ID` = study_id,
    `Index date` = as.character(time_zero),
    `Best RECIST on study` = coalesce(best_recist, "Not assessed"),
    `PFS event` = if_else(pfs_event, "Progression", "Censored"),
    `PFS date` = as.character(pfs_end_date),
    `PFS (months)` = fmt_num(pfs_months, 1),
    `OS event` = if_else(os_event, "Death", "Censored"),
    `OS date` = as.character(os_end_date),
    `OS (months)` = fmt_num(os_months, 1),
    `Final status` = coalesce(final_status, "Not recorded")
  )

ctr_table_style <- function(kbl) {
  kbl %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
    row_spec(0, bold = TRUE, background = "#D9E2F3")
}

# =============================================================================
# R chunk: trial-info-table  ()
# =============================================================================
trial_info <- tibble(
  Field = c(
    "Disease",
    "Stage of disease / treatment",
    "Prior therapy",
    "Type of study",
    "Primary endpoint",
    "Secondary endpoints",
    "Additional details"
  ),
  Description = c(
    "Advanced/recurrent solid tumors (Cohort 1); resected high-risk urothelial carcinoma planned for adjuvant nivolumab (Cohort 2, closed to accrual)",
    "Metastatic/recurrent solid tumors (Cohort 1); post-surgical high-risk urothelial carcinoma (Cohort 2, not enrolled)",
    "Variable prior systemic and local therapy per disease (see Patient Characteristics)",
    "Single-center, open-label, non-randomized pilot",
    "Maximum percentage decline in peripheral blood CRP from cycle 1 day 1 during on-study colchicine",
    "CTCAE v5.0 safety; PFS (Cohort 1); DFS/OS (Cohort 2); exploratory cytokines, ctDNA, tissue correlatives",
    paste0(
      "Cohort 1 low-dose (BID) and high-dose (TID) were sequential; only low-dose enrolled. ",
      "Cohort 2 closed for poor accrual before enrollment. During conduct, the protocol was amended: ",
      "the first ", enrolled_original_n, " enrolled patients received indefinite on-study colchicine ",
      "with serial imaging surveillance; the subsequent ", enrolled_amended_n,
      " enrolled patients received a fixed 2-week colchicine course (phase 0 pharmacodynamic design). ",
      "Exploratory Kaplan–Meier PFS/OS analyses were restricted to the original-schedule cohort (n = ",
      km_eligible_n, ")."
    )
  )
)
trial_info %>%
  kbl(
    caption = tbl_cap("Trial information."),
    col.names = c("Trial information", "Description"),
    align = c("l", "l")
  ) %>%
  ctr_table_style()

# =============================================================================
# R chunk: drug-info-table  ()
# =============================================================================
drug_info <- tibble(
  `Drug information` = c(
    "Generic/working name", "Drug type", "Drug class", "Dose (Cohort 1 low)",
    "Dose (Cohort 1 high)", "Dose (Cohort 2)", "Route", "Schedule"
  ),
  Details = c(
    "Colchicine", "Anti-inflammatory alkaloid", "Microtubule polymerization inhibitor (indirect NLRP3 pathway modulation)",
    "0.6 mg oral BID × 14 days", "0.6 mg oral TID × 14 days (not enrolled)",
    "0.6 mg oral BID × 28 days with nivolumab 480 mg IV every 4 weeks (cohort closed)",
    "Oral", "One cycle per protocol; missed doses not doubled"
  )
)
drug_info %>%
  kbl(
    caption = tbl_cap("Drug information."),
    col.names = c("Drug information", "Details"),
    align = c("l", "l")
  ) %>%
  ctr_table_style()

# =============================================================================
# R chunk: primary-assessment-table  ()
# =============================================================================
primary_assess <- tibble(
  `Primary assessment` = c(
    "Title",
    "Number of patients who signed informed consent",
    "Number of patients enrolled",
    "Number of screen failures before enrollment",
    "Number evaluable for toxicity",
    "Number evaluable for primary efficacy (CRP)",
    "Evaluation method",
    "Primary analysis result (Cohort 1 low dose)"
  ),
  Value = c(
    "Maximum percent decline in CRP from baseline (cycle 1 day 1) during post-baseline C1 on-treatment visits",
    as.character(consented_n),
    as.character(enrolled_n),
    as.character(screen_fail_n),
    as.character(evaluable_toxicity_n),
    as.character(evaluable_efficacy_n),
    "Peripheral blood CRP (mg/L) at protocol C1 visits; baseline = cycle 1 day 1 (or screening fallback if missing); minimum post-baseline value at C1D8/C1D15 used for maximum decline. Cycle 2 values excluded from primary endpoint per protocol treatment window.",
    paste0(
      "Mean max decline ", fmt_num(primary_mean), "% (SD ", fmt_num(primary_sd), "%); ",
      "one-sided one-sample t-test vs 0, t = ", fmt_num(primary_t, 2), ", p = ", fmt_p(primary_p),
      " (prespecified α = 0.025)"
    )
  )
)
primary_assess %>%
  kbl(
    caption = tbl_cap("Primary assessment method."),
    col.names = c("Primary assessment", "Value"),
    align = c("l", "l")
  ) %>%
  ctr_table_style()

# =============================================================================
# R chunk: accrual-disposition-table  ()
# =============================================================================
accrual_overview <- tibble(
  Disposition = c(
    "Signed informed consent",
    "Enrolled on study",
    "  Original schedule (indefinite colchicine + surveillance imaging)",
    "  Amended schedule (fixed 2-week phase 0 intervention)",
    "Screen failure before enrollment",
    "  Before original schedule",
    "  After protocol amendment"
  ),
  N = c(
    consented_n,
    enrolled_n,
    enrolled_original_n,
    enrolled_amended_n,
    screen_fail_n,
    screen_fail_original_n,
    screen_fail_amended_n
  )
) %>%
  mutate(`N (%)` = case_when(
    Disposition == "Signed informed consent" ~ as.character(N),
    Disposition == "Enrolled on study" ~ paste0(N, " (", trimws(fmt_num(100 * N / consented_n, 1)), "% of consented)"),
    Disposition == "Screen failure before enrollment" ~ paste0(N, " (", trimws(fmt_num(100 * N / consented_n, 1)), "% of consented)"),
    grepl("^  Original schedule", Disposition) ~ paste0(N, " (Kaplan–Meier–evaluable)"),
    grepl("^  Amended schedule", Disposition) ~ paste0(N, " (primary CRP / toxicity only)"),
    TRUE ~ paste0(N, " (", trimws(fmt_num(100 * N / consented_n, 1)), "% of consented)")
  ))

accrual_overview %>%
  kbl(
    caption = tbl_cap("Accrual disposition among consented subjects."),
    col.names = c("Disposition", "N", "N (%)"),
    align = c("l", "r", "r")
  ) %>%
  ctr_table_style()

if (nrow(screen_fail_summary) > 0) {
  screen_fail_summary %>%
    mutate(
      `N (%)` = paste0(`N`, " (", trimws(fmt_num(100 * `N` / consented_n, 1)), "% of consented)")
    ) %>%
    select(`Screen-failure reason`, N, `N (%)`, Details) %>%
    kbl(
      caption = tbl_cap("Reasons for screen failure before enrollment."),
      col.names = c("Screen-failure reason", "N", "N (%)", "Details"),
      align = c("l", "r", "r", "l")
    ) %>%
    ctr_table_style()
}

if (screen_fail_n > 0) {
  screen_disposition %>%
    filter(screen_fail) %>%
    transmute(
      `Consent date` = format(consent_date, "%Y-%m-%d"),
      `Screen-failure reason` = screen_fail_reason,
      Details = screen_fail_detail
    ) %>%
    kbl(
      caption = tbl_cap("Screen-failure subjects: reason and supporting detail."),
      col.names = c("Consent date", "Screen-failure reason", "Details"),
      align = c("l", "l", "l")
    ) %>%
    ctr_table_style()
}

# =============================================================================
# R chunk: patient-chars-table  ()
# =============================================================================
n_male <- sum(patient_master$gender_label == "Male", na.rm = TRUE)
n_female <- sum(patient_master$gender_label == "Female", na.rm = TRUE)
n_sex_other <- sum(patient_master$gender_label == "Other", na.rm = TRUE)
n_sex_unknown <- sum(is.na(patient_master$gender_label), na.rm = TRUE)

assert_partition_sums(
  "Patient characteristics: sex",
  c(n_male, n_female, n_sex_other, n_sex_unknown),
  enrolled_n
)

ecog01 <- sum(patient_master$ecog_code %in% c(0, 1), na.rm = TRUE)
ecog2p <- sum(patient_master$ecog_code %in% c(2, 3, 4), na.rm = TRUE)
ecog_unknown <- enrolled_n - ecog01 - ecog2p

assert_partition_sums(
  "Patient characteristics: ECOG",
  c(ecog01, ecog2p, ecog_unknown),
  enrolled_n
)

age_vec <- patient_master$age_years[!is.na(patient_master$age_years)]
age_med <- if (length(age_vec)) median(age_vec) else NA_real_
age_rng <- if (length(age_vec)) range(age_vec) else c(NA_real_, NA_real_)
n_age_known <- length(age_vec)

prior_med <- median(prior_systemic_lines$n_systemic_lines, na.rm = TRUE)
prior_rng <- range(prior_systemic_lines$n_systemic_lines, na.rm = TRUE)

fmt_n_pct <- function(n, d = enrolled_n) paste0(n, " (", trimws(fmt_num(100 * n / d, 1)), "%)")

tumor_type_rows <- tumor_type_summary %>%
  mutate(
    `Patient characteristics` = paste0("  ", primary_tumor_type, ", n (%)"),
    Value = fmt_n_pct(N)
  ) %>%
  select(`Patient characteristics`, Value)

therapy_class_rows <- tibble(therapy_class = therapy_class_levels_table1) %>%
  left_join(
    prior_tx %>%
      group_by(therapy_class) %>%
      summarise(`N patients` = n_distinct(study_id), .groups = "drop"),
    by = "therapy_class"
  ) %>%
  mutate(`N patients` = replace_na(`N patients`, 0L)) %>%
  mutate(
    `Patient characteristics` = paste0("  ", therapy_class, ", n (%)"),
    Value = fmt_n_pct(`N patients`)
  ) %>%
  select(`Patient characteristics`, Value)

pat_chars <- bind_rows(
  tibble(
    `Patient characteristics` = c(
      "Number enrolled (Cohort 1 low dose)",
      "Sex, male, n (%)",
      "Sex, female, n (%)",
      "Age at enrollment, median (range), years",
      "Performance status ECOG 0–1, n (%)",
      "Performance status ECOG ≥2, n (%)",
      "Prior systemic therapy lines, median (range), n"
    ),
    Value = c(
      as.character(enrolled_n),
      fmt_n_pct(n_male),
      fmt_n_pct(n_female),
      if (n_age_known == 0) {
        "Not available"
      } else {
        paste0(
          fmt_num(age_med, 0), " (", fmt_num(age_rng[1], 0), "–", fmt_num(age_rng[2], 0), ")"
        )
      },
      fmt_n_pct(ecog01),
      fmt_n_pct(ecog2p),
      if (length(prior_med) == 0 || is.na(prior_med)) {
        "Not available"
      } else {
        paste0(prior_med, " (", prior_rng[1], "–", prior_rng[2], ")")
      }
    )
  ),
  tibble(`Patient characteristics` = "Exact tumor type", Value = ""),
  tumor_type_rows,
  tibble(`Patient characteristics` = "Prior therapy class (N patients)", Value = ""),
  therapy_class_rows
)
pat_chars %>%
  kbl(
    caption = tbl_cap("Patient characteristics (Cohort 1 low dose)."),
    col.names = c("Characteristic", "Value"),
    align = c("l", "l")
  ) %>%
  ctr_table_style()

# =============================================================================
# R chunk: crp-endpoint-table  ()
# =============================================================================
crp_endpoint %>%
  mutate(
    baseline_crp_mg_l = round(baseline_crp_mg_l, 1),
    min_on_tx_crp_mg_l = round(min_on_tx_crp_mg_l, 1),
    max_pct_decline = round(max_pct_decline, 1)
  ) %>%
  rename(
    `Study ID` = study_id,
    `Baseline CRP (mg/L)` = baseline_crp_mg_l,
    `Lowest on-treatment CRP (mg/L)` = min_on_tx_crp_mg_l,
    `Max % decline` = max_pct_decline,
    `N on-tx CRP values` = n_crp_on_tx
  ) %>%
  kbl(caption = tbl_cap("Primary CRP endpoint derivation by enrolled patient.")) %>%
  ctr_table_style()

# =============================================================================
# R chunk: waterfall-plot  (, fig.cap=fig_cap("Maximum percent decline in C-reactive protein (CRP) from baseline during on-treatment assessment by enrolled patient (Cohort 1 low dose). Abbreviations: CRP, C-reactive protein."))
# =============================================================================
if (nrow(crp_endpoint) > 0) {
  crp_endpoint %>%
    mutate(study_id = reorder(factor(as.character(study_id)), -max_pct_decline)) %>%
    ggplot(aes(x = study_id, y = max_pct_decline)) +
    geom_col(width = 0.75, fill = "#2C6FAC") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      x = "Study ID",
      y = "Maximum CRP decline (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank())
} else {
  plot.new()
  text(0.5, 0.5, "No evaluable CRP endpoint data")
}

# =============================================================================
# R chunk: tumor-waterfall-plot  (, fig.width = 10, fig.height = 7.2, fig.cap=fig_cap("Best percent change in target-lesion sum of diameters (SOD) from baseline among enrolled patients with ≥2 on-study scan sets (n = 5). Axis labels under each Study ID show maximum on-treatment CRP decline (primary endpoint definition; CRP NE if no cycle 1 post-baseline CRP) and the full prior systemic regimen (agents in chronological order). Negative SOD values indicate tumor shrinkage. Abbreviations: SOD, sum of diameters; NE, not evaluable; PLD, pegylated liposomal doxorubicin."))
# =============================================================================
# Full prior systemic regimen (all coded systemic agents, chronological by start date)
systemic_rx_pat <- paste0(
  "(carboplatin|cisplatin|oxaliplatin|taxol|paclitaxel|docetaxel|gemcitabine|",
  "topotecan|doxil|doxorubicin|etoposide|bevacizumab|pembrolizumab|nivolumab|",
  "ipilimumab|atezolizumab|enfortumab|olaparib|rucaparib|niraparib|xeloda|",
  "capecitabine|irinotecan|panitumumab|fluorouracil|5-fu)"
)
normalize_drug <- function(x) {
  dplyr::recode(
    tolower(x),
    taxol = "paclitaxel",
    doxil = "PLD",
    xeloda = "capecitabine",
    "5-fu" = "5-FU",
    fluorouracil = "5-FU",
    .default = tolower(x)
  )
}
prior_full_regimen <- dat %>%
  filter(
    study_id %in% enrolled_ids,
    redcap_repeat_instrument == "prior_treatment_therapy",
    !is.na(prior_treatment_therapies),
    prior_treatment_therapies != ""
  ) %>%
  mutate(
    study_id = as.character(study_id),
    therapy_l = tolower(coalesce(prior_treatment_therapies, "")),
    drug = stringr::str_extract(.data$therapy_l, systemic_rx_pat),
    start_dt = .data$therapy_start_date
  ) %>%
  filter(!is.na(.data$drug)) %>%
  arrange(study_id, .data$start_dt, .data$drug) %>%
  group_by(study_id) %>%
  summarise(
    prior_full = paste(unique(normalize_drug(.data$drug)), collapse = " → "),
    .groups = "drop"
  )

tumor_waterfall_dat <- tumor_change %>%
  mutate(study_id = as.character(study_id)) %>%
  left_join(
    crp_endpoint %>% transmute(study_id = as.character(study_id), max_pct_decline),
    by = "study_id"
  ) %>%
  left_join(prior_full_regimen, by = "study_id") %>%
  mutate(
    plot_y = best_pct_change,
    direction = if_else(best_pct_change < 0, "Shrinkage", "Growth"),
    crp_lab = if_else(
      is.na(max_pct_decline),
      "CRP NE",
      paste0("CRP ↓", sprintf("%.0f", max_pct_decline), "%")
    ),
    prior_lab = if_else(
      is.na(prior_full) | prior_full == "",
      "prior: none coded",
      paste0("prior:\n", stringr::str_wrap(prior_full, width = 22))
    ),
    axis_lab = paste0("ID ", study_id, "\n", crp_lab, "\n", prior_lab)
  ) %>%
  mutate(study_id = reorder(factor(study_id), best_pct_change))

axis_lab_map <- setNames(tumor_waterfall_dat$axis_lab, as.character(tumor_waterfall_dat$study_id))

if (nrow(tumor_waterfall_dat) > 0) {
  ggplot(tumor_waterfall_dat, aes(x = study_id, y = plot_y, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_x_discrete(labels = axis_lab_map) +
    scale_fill_manual(
      values = c("Shrinkage" = "#2C6FAC", "Growth" = "#C75B39"),
      guide = "none"
    ) +
    labs(
      x = NULL,
      y = "Best % change in target-lesion SOD"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 7.5, lineheight = 0.95, vjust = 1),
      plot.margin = margin(5.5, 5.5, 20, 5.5)
    )
} else {
  plot.new()
  text(0.5, 0.5, "No evaluable tumor-change data")
}

# =============================================================================
# R chunk: crp-longitudinal-plot  (, fig.width = 9, fig.height = 11, fig.cap = fig_cap("Longitudinal peripheral blood CRP (mg/L) across protocol C1 visits (baseline C1D1 and post-baseline C1D8/C1D15) for enrolled patients with serial C1 values. Abbreviations: CRP, C-reactive protein."))
# =============================================================================
crp_plot <- crp_long %>%
  mutate(
    visit = case_when(
      redcap_event_name %in% c("cycle_1_day_1__3_d_arm_1", "cycle_1_day_1_arm_2", "cycle_1_day_1_arm_3") ~ "C1D1",
      redcap_event_name %in% c("cycle_1_day_8_arm_2", "cycle_1_day_8_arm_3") ~ "C1D8",
      redcap_event_name %in% c("cycle_1_day_15__3_arm_1", "cycle_1_day_15_arm_2", "cycle_1_day_15_arm_3") ~ "C1D15",
      TRUE ~ NA_character_
    ),
    study_id = factor(study_id),
    visit = factor(visit, levels = c("C1D1", "C1D8", "C1D15"))
  ) %>%
  filter(!is.na(visit)) %>%
  group_by(study_id, visit) %>%
  summarise(crp_value = median(crp_value, na.rm = TRUE), .groups = "drop") %>%
  group_by(study_id) %>%
  filter(n_distinct(visit) >= 2, any(visit %in% c("C1D8", "C1D15"))) %>%
  ungroup()

if (nrow(crp_plot) > 0) {
  ggplot(crp_plot, aes(x = visit, y = crp_value, group = study_id, color = study_id)) +
    geom_line(linewidth = 1.1, alpha = 0.85) +
    geom_point(size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06)), breaks = scales::pretty_breaks(n = 14)) +
    labs(x = "Visit", y = "CRP (mg/L)", color = "Study ID") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")
} else {
  plot.new()
  text(0.5, 0.5, "No CRP longitudinal data")
}

# =============================================================================
# R chunk: toxicity-ae-patient-incidence  (, eval=nrow(ae_by_term_grade_patient) > 0)
# =============================================================================
ae_by_term_grade_patient %>%
  rename(`Adverse event (CTCAE v5)` = ae_term_display, Grade = ae_grade_label) %>%
  mutate(Grade = paste0("Grade ", Grade)) %>%
  select(`Adverse event (CTCAE v5)`, Grade, `N (%)`) %>%
  kbl(caption = tbl_cap(paste0(
    "Treatment-related adverse events by CTCAE v5 preferred term and worst grade per patient (N = ",
    enrolled_n,
    "). Patients with any grade \u22653 event: ",
    round(grade3plus_rate, 0),
    "%."
  ))) %>%
  ctr_table_style()

# =============================================================================
# R chunk: pd-summary  ()
# =============================================================================
pd_long <- dat %>%
  filter(study_id %in% enrolled_ids) %>%
  select(study_id, redcap_event_name, crp_value, ferritin_value_num, esr_value_num) %>%
  pivot_longer(
    cols = c(crp_value, ferritin_value_num, esr_value_num),
    names_to = "marker",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    marker = recode(marker,
      crp_value = "CRP (mg/L)",
      ferritin_value_num = "Ferritin (µg/dL)",
      esr_value_num = "ESR (mm)"
    )
  )

pd_summary <- pd_long %>%
  group_by(marker) %>%
  summarise(
    N = n(),
    Median = median(value),
    IQR_low = quantile(value, 0.25),
    IQR_high = quantile(value, 0.75),
    .groups = "drop"
  ) %>%
  mutate(
    `Median (IQR)` = paste0(
      fmt_num(Median), " (", fmt_num(IQR_low), "–", fmt_num(IQR_high), ")"
    )
  ) %>%
  select(Marker = marker, N, `Median (IQR)`)

pd_summary %>%
  kbl(caption = tbl_cap("Pharmacodynamic laboratory values (all on-study assessments, enrolled patients).")) %>%
  ctr_table_style()

# =============================================================================
# R chunk: pd-change-summary  (, results='asis')
# =============================================================================
# PD change block: baseline (C1D1; screening fallback) to post-baseline C1 (C1D8/C1D15).
pd_marker_cols <- c(
  crp_value = "CRP (mg/L)",
  ferritin_value_num = "Ferritin (µg/dL)",
  esr_value_num = "ESR (mm)"
)

pd_base <- dat %>%
  filter(study_id %in% enrolled_ids, redcap_event_name %in% c(baseline_events, screening_events)) %>%
  select(study_id, redcap_event_name, all_of(names(pd_marker_cols))) %>%
  pivot_longer(
    cols = all_of(names(pd_marker_cols)),
    names_to = "marker_key",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    marker = recode(marker_key, !!!pd_marker_cols),
    source_priority = if_else(redcap_event_name %in% baseline_events, 0L, 1L)
  ) %>%
  arrange(study_id, marker, source_priority, redcap_event_name) %>%
  group_by(study_id, marker) %>%
  slice(1) %>%
  ungroup() %>%
  select(study_id, marker, baseline_value = value)

pd_post <- dat %>%
  filter(study_id %in% enrolled_ids, redcap_event_name %in% post_baseline_events) %>%
  select(study_id, redcap_event_name, all_of(names(pd_marker_cols))) %>%
  pivot_longer(
    cols = all_of(names(pd_marker_cols)),
    names_to = "marker_key",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(marker = recode(marker_key, !!!pd_marker_cols)) %>%
  group_by(study_id, marker) %>%
  summarise(
    post_min_value = min(value, na.rm = TRUE),
    n_post_values = n(),
    .groups = "drop"
  )

pd_change <- pd_base %>%
  inner_join(pd_post, by = c("study_id", "marker")) %>%
  filter(!is.na(baseline_value), baseline_value > 0) %>%
  mutate(
    pct_change = 100 * (post_min_value - baseline_value) / baseline_value,
    pct_decline = pmax(-pct_change, 0)
  )

pd_change_summary <- pd_change %>%
  group_by(marker) %>%
  summarise(
    `N evaluable` = n_distinct(study_id),
    `Baseline median` = fmt_num(median(baseline_value, na.rm = TRUE), 1),
    `Post-baseline minimum median` = fmt_num(median(post_min_value, na.rm = TRUE), 1),
    `Median % decline (IQR)` = paste0(
      fmt_num(median(pct_decline, na.rm = TRUE), 1), " (",
      fmt_num(quantile(pct_decline, 0.25, na.rm = TRUE), 1), "–",
      fmt_num(quantile(pct_decline, 0.75, na.rm = TRUE), 1), ")"
    ),
    .groups = "drop"
  ) %>%
  rename(Marker = marker)

pd_change_summary %>%
  kbl(caption = tbl_cap("Pharmacodynamic decline from baseline to minimum post-baseline C1 value (C1D8/C1D15). Decline is bounded at 0% (no decline).")) %>%
  ctr_table_style()

ferritin_c1_evaluable_n <- pd_change %>%
  filter(marker == "Ferritin (µg/dL)") %>%
  summarise(n = n_distinct(study_id)) %>%
  pull(n)

ferritin_any_post_evaluable_n <- {
  pd_post_any <- dat %>%
    filter(study_id %in% enrolled_ids) %>%
    filter(!redcap_event_name %in% c(screening_events, baseline_events)) %>%
    select(study_id, redcap_event_name, ferritin_value_num) %>%
    filter(!is.na(ferritin_value_num)) %>%
    distinct(study_id)

  pd_base_ferritin <- pd_base %>%
    filter(marker == "Ferritin (µg/dL)") %>%
    distinct(study_id)

  inner_join(pd_base_ferritin, pd_post_any, by = "study_id") %>%
    summarise(n = n_distinct(study_id)) %>%
    pull(n)
}

cat(
  "\n\n**Ferritin evaluable count note:** baseline-to-post C1 (D8/D15) ferritin was available for ",
  ferritin_c1_evaluable_n,
  " patient(s); using any post-baseline visit (any cycle) would yield ",
  ferritin_any_post_evaluable_n,
  " patient(s).\n\n",
  sep = ""
)

# =============================================================================
# R chunk: pd-crp-responders  ()
# =============================================================================
crp_response <- pd_change %>%
  filter(marker == "CRP (mg/L)") %>%
  summarise(
    `N evaluable` = n(),
    `Decline ≥30%, n (%)` = paste0(sum(pct_decline >= 30), " (", sprintf("%.1f", 100 * mean(pct_decline >= 30)), "%)"),
    `Decline ≥50%, n (%)` = paste0(sum(pct_decline >= 50), " (", sprintf("%.1f", 100 * mean(pct_decline >= 50)), "%)"),
    `Median CRP % decline (IQR)` = paste0(
      fmt_num(median(pct_decline, na.rm = TRUE), 1), " (",
      fmt_num(quantile(pct_decline, 0.25, na.rm = TRUE), 1), "–",
      fmt_num(quantile(pct_decline, 0.75, na.rm = TRUE), 1), ")"
    )
  )

crp_response %>%
  kbl(caption = tbl_cap("CRP pharmacodynamic response summary among evaluable patients.")) %>%
  ctr_table_style()

# =============================================================================
# R chunk: bor-summary  ()
# =============================================================================
bor_tbl <- surv_dat %>%
  transmute(`Best overall response (RECIST 1.1)` = coalesce(best_recist, "Not assessed")) %>%
  count(`Best overall response (RECIST 1.1)`, name = "N patients") %>%
  mutate(`N (%)` = paste0(`N patients`, " (", sprintf("%.1f", 100 * `N patients` / km_eligible_n), "%)")) %>%
  arrange(desc(`N patients`), `Best overall response (RECIST 1.1)`)

bor_tbl %>%
  kbl(caption = tbl_cap("Exploratory best overall response (RECIST 1.1) among original-schedule enrolled patients (Kaplan–Meier cohort).")) %>%
  ctr_table_style()

# =============================================================================
# R chunk: pfs-os-summary  ()
# =============================================================================
surv_summary <- tibble(
  Endpoint = c("Progression-free survival", "Overall survival"),
  `N patients` = c(nrow(surv_pfs), nrow(surv_os)),
  `N events` = c(sum(surv_pfs$pfs_event), sum(surv_os$os_event)),
  `Median (months, 95% CI)` = c(
    if (is.na(pfs_med$median)) {
      "Not reached / not estimable"
    } else {
      paste0(
        fmt_num(pfs_med$median, 1), " (",
        fmt_num(pfs_med$lcl, 1), "–", fmt_num(pfs_med$ucl, 1), ")"
      )
    },
    if (is.na(os_med$median)) {
      "Not reached / not estimable"
    } else {
      paste0(
        fmt_num(os_med$median, 1), " (",
        fmt_num(os_med$lcl, 1), "–", fmt_num(os_med$ucl, 1), ")"
      )
    }
  )
)
surv_summary %>%
  kbl(caption = tbl_cap(paste0(
    "Exploratory progression-free survival and overall survival summary (Kaplan–Meier medians; original-schedule cohort, n = ",
    km_eligible_n, ")."
  ))) %>%
  ctr_table_style()

# =============================================================================
# R chunk: pfs-os-by-patient  ()
# =============================================================================
pfs_os_tbl %>%
  kbl(
    caption = tbl_cap("Patient-level progression-free survival and overall survival among original-schedule enrolled patients (derived from RECIST, discontinuation, and final status fields).")
  ) %>%
  ctr_table_style()

# =============================================================================
# R chunk: km-pfs  (, fig.cap=fig_cap(paste0("Progression-free survival (Kaplan–Meier) from treatment start among patients enrolled under the original indefinite colchicine schedule with surveillance imaging (n = ", km_eligible_n, "; exploratory; not powered). PFS events: RECIST PD or discontinuation for progression. Abbreviations: PFS, progression-free survival; PD, progressive disease; RECIST, Response Evaluation Criteria in Solid Tumors.")), eval=!is.null(fit_pfs))
# =============================================================================
if (!is.null(fit_pfs)) {
  ggsurvplot(
    fit_pfs,
    data = surv_pfs,
    risk.table = TRUE,
    conf.int = TRUE,
    xlab = "Months from treatment start",
    ylab = "Progression-free survival probability",
    title = "Exploratory progression-free survival",
    legend = "none",
    risk.table.height = 0.25
  )
}

# =============================================================================
# R chunk: km-os  (, fig.cap=fig_cap(paste0("Overall survival (Kaplan–Meier) from treatment start among patients enrolled under the original indefinite colchicine schedule with surveillance imaging (n = ", km_eligible_n, "; exploratory; not powered). Abbreviations: OS, overall survival.")), eval=!is.null(fit_os))
# =============================================================================
if (!is.null(fit_os)) {
  ggsurvplot(
    fit_os,
    data = surv_os,
    risk.table = TRUE,
    conf.int = TRUE,
    xlab = "Months from treatment start",
    ylab = "Overall survival probability",
    title = "Exploratory overall survival",
    legend = "none",
    risk.table.height = 0.25
  )
}

# =============================================================================
# R chunk: appendix-included-subjects  ()
# =============================================================================
included_subject_ids <- sort(unique(crp_endpoint$study_id))

included_crp_visits <- crp_long %>%
  filter(study_id %in% included_subject_ids) %>%
  mutate(
    visit = case_when(
      redcap_event_name %in% c("cycle_1_day_1__3_d_arm_1", "cycle_1_day_1_arm_2", "cycle_1_day_1_arm_3") ~ "C1D1",
      redcap_event_name %in% c("cycle_1_day_8_arm_2", "cycle_1_day_8_arm_3") ~ "C1D8",
      redcap_event_name %in% c("cycle_1_day_15__3_arm_1", "cycle_1_day_15_arm_2", "cycle_1_day_15_arm_3") ~ "C1D15",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(visit)) %>%
  group_by(study_id, visit) %>%
  summarise(crp_mg_l = min(crp_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = visit, values_from = crp_mg_l, names_prefix = "CRP_")

included_tbl <- tibble(study_id = included_subject_ids) %>%
  left_join(
    patient_master %>%
      select(
        study_id,
        primary_tumor_type,
        age_years,
        gender_label,
        ecog_label
      ),
    by = "study_id"
  ) %>%
  left_join(prior_systemic_lines, by = "study_id") %>%
  left_join(
    crp_endpoint %>%
      select(study_id, baseline_crp_mg_l, min_on_tx_crp_mg_l, max_pct_decline, n_crp_on_tx),
    by = "study_id"
  ) %>%
  left_join(included_crp_visits, by = "study_id") %>%
  left_join(
    surv_dat %>%
      select(study_id, pfs_event, pfs_months, os_event, os_months),
    by = "study_id"
  ) %>%
  mutate(
    `Study ID` = study_id,
    `Primary cancer type` = primary_tumor_type,
    `Age (years)` = round(age_years, 0),
    Sex = gender_label,
    ECOG = ecog_label,
    `Prior systemic lines` = n_systemic_lines,
    `CRP baseline C1D1 (mg/L)` = round(baseline_crp_mg_l, 1),
    `CRP C1D8 (mg/L)` = round(CRP_C1D8, 1),
    `CRP C1D15 (mg/L)` = round(CRP_C1D15, 1),
    `Lowest post-baseline CRP (mg/L)` = round(min_on_tx_crp_mg_l, 1),
    `Max CRP decline (%)` = round(max_pct_decline, 1),
    `N post-baseline CRP values` = n_crp_on_tx,
    `PFS event` = if_else(pfs_event, "Yes", "No"),
    `PFS (months)` = round(pfs_months, 1),
    `OS event (death)` = if_else(os_event, "Yes", "No"),
    `OS (months)` = round(os_months, 1)
  ) %>%
  select(
    `Study ID`,
    `Primary cancer type`,
    `Age (years)`,
    Sex,
    ECOG,
    `Prior systemic lines`,
    `CRP baseline C1D1 (mg/L)`,
    `CRP C1D8 (mg/L)`,
    `CRP C1D15 (mg/L)`,
    `Lowest post-baseline CRP (mg/L)`,
    `Max CRP decline (%)`,
    `N post-baseline CRP values`,
    `PFS event`,
    `PFS (months)`,
    `OS event (death)`,
    `OS (months)`
  ) %>%
  arrange(`Study ID`)

included_tbl %>%
  kbl(
    caption = tbl_cap(paste0(
      "Subject-level dataset for patients included in primary endpoint analysis (N = ",
      nrow(included_tbl),
      ")."
    ))
  ) %>%
  ctr_table_style()
