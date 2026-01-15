library(dplyr)
library(clubSandwich)

measure_group <- function(measure, effect_is_log) {
  if (is.character(effect_is_log)) {
    effect_is_log <- tolower(effect_is_log) == "true"
  }
  if (is.factor(effect_is_log)) {
    effect_is_log <- tolower(as.character(effect_is_log)) == "true"
  }
  if (!is.na(measure)) {
    if (measure %in% c("RR", "OR", "HR", "PETO")) return("log_ratio")
    if (measure %in% c("MD", "SMD", "RD")) return("mean_diff")
    if (measure %in% c("GEN")) return("generic")
  }
  if (!is.na(effect_is_log) && effect_is_log) return("log_ratio")
  "other"
}

categorize_analysis_name <- function(name) {
  if (is.na(name) || is.null(name)) return("other")
  x <- tolower(name)
  if (grepl("mortality|death|survival|fatal", x)) return("mortality")
  if (grepl("adverse|toxicity|side effect|complication|bleed|bleeding|haemorrhage|hemorrhage|infection|sepsis|harm|safety", x)) return("adverse_events")
  if (grepl("pain|analges", x)) return("pain")
  if (grepl("quality of life|qol|well[- ]?being|hrqol", x)) return("quality_of_life")
  if (grepl("length of stay|hospital stay|days in hospital|hospitalization", x)) return("length_of_stay")
  if (grepl("relapse|recurrence|progression", x)) return("relapse")
  if (grepl("response|remission|cure|recovery|resolution", x)) return("response")
  if (grepl("score|index|scale|rating|questionnaire|functional|disability", x)) return("function_score")
  if (grepl("blood pressure|hb|hba1c|cholesterol|lipid|glucose|viral load|titer|titre|concentration|level", x)) return("biomarker")
  if (grepl("cost|economic|resource|budget", x)) return("cost")
  "other"
}

input_path <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output/sensitivity_fixed_influence.csv"
output_path <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output/advanced_mixed_effects_logit.csv"

df <- read.csv(input_path)
analysis_df <- df %>%
  filter(k_post >= 2, k_pre > 0, !is.na(delta_i2_re)) %>%
  mutate(
    outcome_category = vapply(analysis_name, categorize_analysis_name, character(1)),
    measure_group = mapply(measure_group, measure, effect_is_log),
    log_k_all = log1p(k_all)
  )

category_counts <- table(analysis_df$outcome_category)
small_cats <- names(category_counts[category_counts < 30])
analysis_df$outcome_category <- ifelse(analysis_df$outcome_category %in% small_cats, "other_small", analysis_df$outcome_category)
analysis_df$delta_i2_gt0 <- analysis_df$delta_i2_re > 0

fit_glm <- glm(delta_i2_gt0 ~ frac_pre + log_k_all + data_type + measure_group + outcome_category,
              data = analysis_df, family = binomial())
vc_g <- vcovCR(fit_glm, cluster = analysis_df$dataset_id, type = "CR2")
coef_tab_g <- coef_test(fit_glm, vcov = vc_g, test = "Satterthwaite")
write.csv(as.data.frame(coef_tab_g), output_path, row.names = FALSE)
