# STUDY-21-001175 — prespecified and exploratory analyses (R 4.x)
# Data cut: ProtocolNo2101175PRM_DATA_2026-05-18_2135.csv
# Cohort 1 low dose only (colchicine 0.6 mg oral BID x 14 days)

library(dplyr)
library(survival)

dat <- read.csv("ProtocolNo2101175PRM_DATA_2026-05-18_2135.csv")

# --- Primary endpoint: max % CRP decline from C1D1 during C1 on-treatment visits ---
crp_c1 <- dat %>%
  filter(redcap_event_name %in% c("cycle_1_day_1_arm_1", "cycle_1_day_8_arm_1", "cycle_1_day_15_arm_1")) %>%
  select(study_id, redcap_event_name, c_reactive_protein_mg_l) %>%
  tidyr::pivot_wider(names_from = redcap_event_name, values_from = c_reactive_protein_mg_l)

crp_eval <- crp_c1 %>%
  mutate(
    baseline_crp = cycle_1_day_1_arm_1,
    post_min_crp = pmin(cycle_1_day_8_arm_1, cycle_1_day_15_arm_1, na.rm = TRUE),
    max_crp_decline_pct = 100 * (baseline_crp - post_min_crp) / baseline_crp
  ) %>%
  filter(!is.na(max_crp_decline_pct), is.finite(max_crp_decline_pct))

# Prespecified one-sided one-sample t-test (H0: mean decline = 0), alpha = 0.025
primary_test <- t.test(crp_eval$max_crp_decline_pct, mu = 0, alternative = "greater")
primary_test  # mean 22.6%, SD 32.9%, t = 1.94, p = 0.047 (N = 8 evaluable)

# --- Exploratory Kaplan-Meier PFS/OS (not powered) ---
# PFS event: RECIST PD or discontinuation for progression
# OS event: death from any cause

pfs_os <- dat %>%
  group_by(study_id) %>%
  summarise(
    pfs_time_months = first(pfs_months),
    pfs_event = first(as.integer(pfs_event_flag)),
    os_time_months = first(os_months),
    os_event = first(as.integer(os_event_flag)),
    .groups = "drop"
  )

fit_pfs <- survfit(Surv(pfs_time_months, pfs_event) ~ 1, data = pfs_os)
fit_os  <- survfit(Surv(os_time_months, os_event) ~ 1, data = pfs_os)
summary(fit_pfs)
summary(fit_os)

# Post-hoc power for observed primary effect (N = 8, alpha = 0.025 one-sided)
pwr::power.t.test(n = 8, delta = 22.6 / 32.9, sd = 1, sig.level = 0.025, type = "one.sample", alternative = "one.sided")
