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
#'
#' Also returns CRP at the imaging timepoint corresponding to the best SOD change:
#' post-baseline scan sets are matched in order to sorted target-lesion/BOR assessment
#' dates, and CRP is taken as the closest CRP date within \code{crp_window_days}.
build_protocol2101175_tumor_change <- function(dat, enrolled_ids, crp_window_days = 21) {
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

  empty <- tibble::tibble(
    study_id = character(),
    baseline_sod = numeric(),
    best_pct_change = numeric(),
    max_shrinkage_pct = numeric(),
    n_post_scans = integer(),
    best_scan_set = integer(),
    imaging_date = as.Date(character()),
    baseline_crp = numeric(),
    crp_at_imaging = numeric(),
    crp_date_at_imaging = as.Date(character()),
    crp_pct_change_at_imaging = numeric()
  )

  if (nrow(les) == 0) return(empty)

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

  sod_long <- les %>%
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
    dplyr::ungroup()

  tumor <- sod_long %>%
    dplyr::filter(!is.na(.data$pct_change)) %>%
    dplyr::group_by(.data$study_id) %>%
    dplyr::summarise(
      baseline_sod = dplyr::first(.data$baseline_sod),
      best_pct_change = min(.data$pct_change, na.rm = TRUE),
      max_shrinkage_pct = pmax(-min(.data$pct_change, na.rm = TRUE), 0),
      n_post_scans = dplyr::n(),
      best_scan_set = .data$scan_set[which.min(.data$pct_change)][[1]],
      .groups = "drop"
    )

  if (nrow(tumor) == 0) return(empty)

  # Assessment dates: target-lesion eval and BOR dates, sorted per patient
  assess_dates <- dat %>%
    dplyr::filter(.data$study_id %in% tumor$study_id) %>%
    dplyr::mutate(study_id = as.character(.data$study_id)) %>%
    dplyr::transmute(
      study_id,
      assess_date = dplyr::coalesce(
        .data$date_of_target_lesions_eva,
        .data$date_of_overall_response
      )
    ) %>%
    dplyr::filter(!is.na(.data$assess_date)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$study_id, .data$assess_date) %>%
    dplyr::group_by(.data$study_id) %>%
    dplyr::mutate(post_scan_idx = dplyr::row_number()) %>%
    dplyr::ungroup()

  # Map best post-baseline scan_set (2,3,...) -> 1st, 2nd, ... assessment date
  tumor <- tumor %>%
    dplyr::mutate(post_scan_idx = as.integer(.data$best_scan_set) - 1L) %>%
    dplyr::left_join(assess_dates, by = c("study_id", "post_scan_idx")) %>%
    dplyr::rename(imaging_date = .data$assess_date)

  baseline_events <- c(
    "cycle_1_day_1__3_d_arm_1", "cycle_1_day_1_arm_2", "cycle_1_day_1_arm_3"
  )
  screening_events <- c(
    "screening_28_days_arm_1", "screening_arm_2", "screening_arm_3"
  )

  crp_rows <- dat %>%
    dplyr::filter(.data$study_id %in% tumor$study_id) %>%
    dplyr::mutate(
      study_id = as.character(.data$study_id),
      crp_value = suppressWarnings(as.numeric(.data$c_reactive_protein_value)),
      crp_date = .data$date_c_reactive_protein_co
    ) %>%
    dplyr::filter(!is.na(.data$crp_value), !is.na(.data$crp_date)) %>%
    dplyr::select(study_id, redcap_event_name, crp_value, crp_date)

  pick_baseline_crp <- function(id) {
    sub <- crp_rows %>% dplyr::filter(.data$study_id == id)
    base <- sub %>%
      dplyr::filter(.data$redcap_event_name %in% baseline_events) %>%
      dplyr::arrange(.data$crp_date) %>%
      dplyr::slice(1)
    if (nrow(base) == 0) {
      base <- sub %>%
        dplyr::filter(.data$redcap_event_name %in% screening_events) %>%
        dplyr::arrange(.data$crp_date) %>%
        dplyr::slice(1)
    }
    if (nrow(base) == 0) {
      return(list(value = NA_real_, date = as.Date(NA)))
    }
    list(value = base$crp_value[[1]], date = base$crp_date[[1]])
  }

  pick_crp_at_imaging <- function(id, imaging_date) {
    if (is.na(imaging_date)) {
      return(list(value = NA_real_, date = as.Date(NA)))
    }
    sub <- crp_rows %>% dplyr::filter(.data$study_id == id)
    if (nrow(sub) == 0) {
      return(list(value = NA_real_, date = as.Date(NA)))
    }
    sub <- sub %>%
      dplyr::mutate(abs_days = abs(as.numeric(.data$crp_date - imaging_date))) %>%
      dplyr::filter(.data$abs_days <= crp_window_days) %>%
      dplyr::arrange(.data$abs_days, .data$crp_date)
    if (nrow(sub) == 0) {
      return(list(value = NA_real_, date = as.Date(NA)))
    }
    list(value = sub$crp_value[[1]], date = sub$crp_date[[1]])
  }

  tumor$baseline_crp <- NA_real_
  tumor$crp_at_imaging <- NA_real_
  tumor$crp_date_at_imaging <- as.Date(NA)
  for (i in seq_len(nrow(tumor))) {
    b <- pick_baseline_crp(tumor$study_id[[i]])
    tumor$baseline_crp[[i]] <- b$value
    cimg <- pick_crp_at_imaging(tumor$study_id[[i]], tumor$imaging_date[[i]])
    tumor$crp_at_imaging[[i]] <- cimg$value
    tumor$crp_date_at_imaging[[i]] <- cimg$date
  }

  tumor %>%
    dplyr::mutate(
      crp_pct_change_at_imaging = dplyr::if_else(
        !is.na(.data$baseline_crp) & .data$baseline_crp > 0 & !is.na(.data$crp_at_imaging),
        100 * (.data$baseline_crp - .data$crp_at_imaging) / .data$baseline_crp,
        NA_real_
      ),
      study_id = as.character(.data$study_id)
    ) %>%
    dplyr::select(
      study_id, baseline_sod, best_pct_change, max_shrinkage_pct, n_post_scans,
      best_scan_set, imaging_date, baseline_crp, crp_at_imaging, crp_date_at_imaging,
      crp_pct_change_at_imaging
    ) %>%
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
