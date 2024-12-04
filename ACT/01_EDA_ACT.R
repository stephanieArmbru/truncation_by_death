################################### ACT data ###################################
########################### Explanatory Data Analysis ##########################

##### Stephanie Armbruster
##### August 2024

# Data exploration
# number of deaths among AD and non-AD patients
# number of AD diagnosis
# time of AD diagnosis / death

# CASI score development, conditional on AD diagnosis: barplot over time
# (or age) colored according to AD diagnosis and truncated by death


library("devtools")
install_github("harrisonreeder/SemiCompRisksFreq")

# Libraries ---------------------------------------------------------------
library(haven)
library(SemiCompRisksFreq)
library(survival)
library(tidyverse)


# Data --------------------------------------------------------------------
#read in full visit-level data
act_visit <- haven::read_sas(data_file = "/Users/stephaniearmbruster/Desktop/Grant/Truncation\ by\ Death/ACT/Input/bienvisit.sas7bdat")
# save original copy
act_visit_old <- act_visit

#look at all the variables
summary(act_visit)

# first visit
min(act_visit$visitdt)
min(act_visit$deathdt,na.rm=TRUE)

#it appears the data were pulled on or around 2020-03-05, so we can use that as the
#'administrative censoring' date
max(act_visit$visitdt)
max(act_visit$deathdt,na.rm=TRUE)

# Demographics ------------------------------------------------------------
# censoring post age 99
act_visit <- act_visit %>%
  filter(age < 100)

# set age at death NA if > 99
act_visit %>%
  filter(Age_Death > 99, VISIT == max(VISIT)) %>%
  nrow()

#make table just of baseline visits to have person-level info
act_person <- act_visit %>%
  group_by(SUBJECT) %>%
  filter(VISIT==0)

act_visit <- act_visit %>%
  mutate(Age_Death = ifelse(Age_Death > 99, NA, Age_Death))
act_visit$Age_Death %>% na.omit() %>% range()
act_visit$deathdt %>% na.omit() %>%  range()

# number patients
N <- act_visit$SUBJECT %>% unique() %>% length()
# x[!(x %in% unique(act_person$SUBJECT))]
# act_person %>% filter(!(x %in% SUBJECT))

# gender
act_person$gender %>% table() %>% prop.table()

# average age at baseline
act_person$age %>% mean()
act_visit$age %>% range()


# race
act_person$race %>% table() %>% prop.table()

# average years of visits since baseline
act_visit$VISIT %>% table()

# number of visits
numb_visits <- act_visit %>%
  mutate(anydementia = ifelse(is.na(anydementia), 0, anydementia),
         anyad = ifelse(is.na(anyad), 0, anyad)) %>%
  group_by(SUBJECT) %>%
  summarize(number_visits = length(visitdt),
            recruitment_day = min(visitdt),
            death_date = min(deathdt),
            has_Dementia = max(anydementia),
            had_AD = max(anyad))

numb_visits$number_visits %>% summary()
numb_visits$number_visits %>% range()
numb_visits$number_visits %>% mean()

act_person$deathdt %>% is.na() %>% table() %>% prop.table()
act_person$Age_Death %>% is.na() %>% table()

numb_visits %>% filter(number_visits == 1) %>% nrow()
numb_visits %>% filter(number_visits == 1) %>% pull(recruitment_day) %>% summary()
numb_visits %>% pull(recruitment_day) %>% summary()

# APOE allel
act_person %>% pull(apoe) %>% table(useNA = "always")
# %>% prop.table()

# smoking
# group current and former smokers together
act_visit <- act_visit %>%
  mutate(smoker = case_when(smoke == 0 ~ "never smoker",
                            is.na(smoke) ~ "no info",
                            .default = "former / current smoker"))

act_person$smoke %>% table(useNA = "always")

# reporting smoke status unclear for some subjects
act_visit %>%
  group_by(SUBJECT) %>%
  summarize(smoker = paste0(unique(smoker), collapse = ", ")) %>%
  pull(smoker) %>%
  table()

# anyone who ever said they smoked is considered a former / current smoker
act_visit %>%
  group_by(SUBJECT) %>%
  summarize(smoker = any(smoker == "former / current smoker")) %>%
  ungroup() %>%
  pull(smoker) %>%
  table(useNA = "always")

act_visit %>%
  group_by(SUBJECT) %>%
  summarize(never_smoker = all(smoker == "never smoker" | smoker == "no info") & any(smoker == "never smoker")) %>%
  ungroup() %>%
  pull(never_smoker) %>%
  table(useNA = "always")

# adapt smoking indicator accordingly
act_visit <- act_visit %>%
  group_by(SUBJECT) %>%
  mutate(smoker_ind = case_when(any(smoker == "former / current smoker") ~ "former / current smoker",
                                never_smoker = all(smoker == "never smoker" | smoker == "no info") & any(smoker == "never smoker") ~ "never smoker",
                                .default = NA))

act_visit %>%
  filter(VISIT == max(VISIT)) %>%
  pull(smoker_ind) %>%
  table(useNA = "always")


# Administrative censoring ------------------------------------------------
numb_visits %>% filter(number_visits == 1) %>% nrow()
numb_visits %>% filter(number_visits == 1) %>% pull(recruitment_day) %>% summary()
numb_visits %>% pull(recruitment_day) %>% summary()
# people with one visit seem to have joined later (administrative censoring)

numb_visits %>% filter(number_visits == 1) %>% pull(death_date) %>% is.na() %>% table()
numb_visits %>% filter(number_visits == 1) %>% pull(has_Dementia) %>% table()
numb_visits %>% filter(number_visits == 1) %>% pull(had_AD) %>% table()


# CASI score --------------------------------------------------------------
# 26,056 valid CASI scores among all participants
table(act_visit$casi_valid, useNA = "always")

table(act_visit %>%
        filter(CASI_SC >= 101 | is.na(CASI_SC)) %>%
        pull(CASI_SC),
      useNA = "always")
# 520 are invalid, but only 506 have invalid scores
(act_visit$casi_valid != 1) %>% na.omit() %>% sum()
(act_visit$CASI_SC >= 101) %>% na.omit() %>% sum()

# minimum is 14, maximum is 100 (threshold 85)
act_visit %>% filter(CASI_SC <= 100) %>% pull(CASI_SC) %>% range()

# correct for invalidity (CASI_SC 1000 for invalid according to casi_valid)
act_visit <- act_visit %>%
  mutate(CASI_SC = ifelse(casi_valid != 1 & CASI_SC <= 100 & !is.na(CASI_SC),
                          1000, CASI_SC))

# 5 are valid but have invalid CASI score
act_visit %>%
  filter(CASI_SC > 100 & casi_valid == 1) %>%
  dplyr::select(CASI_SC, casi_valid)

# correct for invalidity (casi_valid 10 for invalid according to CASI_SC)
act_visit <- act_visit %>%
  mutate(casi_valid = ifelse(CASI_SC > 100 & casi_valid == 1,
                             10,
                             casi_valid))

# 525 invalid CASI scores
act_visit %>% filter(casi_valid != 1 & CASI_SC > 100) %>% nrow()
act_visit %>% filter(casi_valid != 1 | CASI_SC > 100) %>% nrow()

# binary indicator for validity of CASI score
# continuous variable for valid CASI score; NA for invalid
act_visit <- act_visit %>%
  mutate(casi_is_valid = ifelse(casi_valid == 1 & CASI_SC <= 100,
                                TRUE,
                                FALSE)) %>%
  mutate(casi_score = ifelse(casi_is_valid == 1,
                             CASI_SC,
                             NA))

act_visit$casi_is_valid %>% table()



# Death -------------------------------------------------------------------
#check to make sure there were no visits that happened after death (there is one!!)
# View(act_visit[which(act_visit$visitdt > act_visit$deathdt),])
# View(act_visit %>% filter(SUBJECT == 154130))

act_visit %>% dplyr::select(deathdt, SUBJECT)
act_visit %>% dplyr::select(Age_Death, SUBJECT)

act_visit %>% nrow()
act_visit %>% filter(visitdt < deathdt) %>% nrow()
act_visit %>% filter(visitdt == deathdt) %>% nrow()

# binary indicator for death process
act_visit <- act_visit %>%
  ungroup() %>%
  mutate(c4dth = ifelse(visitdt < deathdt | is.na(deathdt),
                        0,
                        1)) %>%
  group_by(SUBJECT) %>%
  mutate(c4dth = ifelse(VISIT == max(VISIT) & !is.na(deathdt),
                        1,
                        c4dth))

act_visit$c4dth %>% table()
3117 / N

# binary indicator for dying at some point in the observation period
act_visit <- act_visit %>%
  mutate(dies = !is.na(Age_Death))

(!is.na(act_person$deathdt)) %>% sum()
(!is.na(act_person$Age_Death)) %>% sum()

# Kaplan Meier estimate for survival
act_KM <- act_visit %>%
  group_by(SUBJECT) %>%
  filter(row_number() == max(row_number())) %>%
  mutate(t2death = ifelse(!is.na(Age_Death), Age_Death, age),
         c4death = ifelse(!is.na(Age_Death), 1, 0))


KM_est <- survfit(Surv(time = act_KM$t2death,
                       event = act_KM$c4death) ~ 1)

ggplot(mapping = aes(x = KM_est$time,
           y = KM_est$surv)) +
  geom_line() +
  theme_bw() +
  labs(x = "age",
       y = "Survival probability")



# Dementia / AD  ----------------------------------------------------------
# missing values if test not performed (CASI score < 85) or died or no MD visit
with(act_visit, table(anydementia|anyad, useNA = "always"))

#look at ad/dementia diagnoses by visit
with(act_visit, table(VISIT, anydementia|anyad , useNA = "always"))


#look at correspondence of dementia flag and dsmivdx diagnostic categorical
with(act_visit, table(dsmivdx, anydementia, useNA = "always"))

#make sure everyone with a dementia onset date falls into a dementia category
with(act_visit, table(is.na(onsetdate), dsmivdx, useNA = "always"))
with(act_visit, table(is.na(onsetdate), anyad, useNA = "always"))
with(act_visit, table(is.na(onsetdate), anydementia, useNA = "always"))
#make sure everyone with a dementia age falls into a dementia category
with(act_visit, table(is.na(Age_Onset), dsmivdx, useNA = "always"))


#look at correspondence of ad flag and nindx diagnostic categorical
with(act_visit, table(nindx, anyad, useNA = "always"))

#look at correspondence of dementia and nindx diagnostic categoricals (not complete overlap!)
with(act_visit, table(nindx, dsmivdx, useNA = "always"))
# 23 observations where no dementia but AD diagnosis, those are missing onset dates
with(act_visit, table(anydementia, anyad, useNA = "always"))
# 1 observation with DSMIV diagnoses no dementia, but NINCDS diagnoses dementia
act_visit %>%
  filter(dsmivdx == 0 & nindx == 3) %>%
  dplyr::select(anydementia, anyad)


# date of AD / dementia onset
# act_visit %>% dplyr::select(onsetdate, deathdt, visitdt, SUBJECT) %>% view()
# act_visit %>% dplyr::select(Age_Onset, SUBJECT) %>% view()

# no further observations after dementia onset
# truncation by death and by dementia onset
d1 <- act_visit %>%
  group_by(SUBJECT) %>%
  filter(visitdt == max(visitdt) & anydementia == 1)
d1 %>% nrow()

d2 <- act_visit %>%
  group_by(SUBJECT) %>%
  filter(!is.na(onsetdate))
d2 %>% nrow()

identical(d1, d2)


# counting process indicator for dementia and AD
act_visit <- act_visit %>%
  group_by(SUBJECT) %>%
  mutate(c4Dementia = ifelse(!is.na(onsetdate),
                             1,
                             0))
# counting process indicator for AD
act_visit <- act_visit %>%
  group_by(SUBJECT) %>%
  mutate(c4AD = ifelse(anydementia == 1 & anyad == 1 & !is.na(onsetdate),
                       1,
                       0))


# date of earliest dementia diagnosis
act_visit$onsetdate %>% na.omit() %>% range()

# number of dementia diagnosis
act_visit$anydementia %>% table(useNA = "always")
act_visit$onsetdate %>% is.na() %>% table()
1337 / nrow(act_person)

# average age of dementia diagnosis
act_visit %>% filter(c4Dementia == 1) %>% pull(Age_Onset) %>% mean()

# number of AD diagnosis
act_visit$c4AD %>% table()
1092 / nrow(act_person)

# average age of AD diagnosis
act_visit %>% filter(c4AD == 1) %>% pull(Age_Onset) %>% mean()

# create competing risk
act_visit <- act_visit %>%
  mutate(c4ADDementia = case_when(c4AD == 1 & dies == 1 ~ "AD - dead",
                                  c4AD == 1 & dies == 0 ~ "AD - alive",
                                  c4Dementia == 1 & dies == 1 ~ "other dementia - dead",
                                  c4Dementia == 1 & dies == 0 ~ "other dementia - alive",
                                  c4AD == 0 & c4Dementia == 0 & dies == 1 ~ "no dementia - dead",
                                  c4AD == 0 & c4Dementia == 0 & dies == 0 ~ "no dementia - alive"
                                  )) %>%
  mutate(ADDementia_color = case_when(c4AD == 1  ~ "AD",
                                      c4Dementia == 1 ~ "other dementia",
                                      c4AD == 0 & c4Dementia == 0 & dies == 1 ~ "Truncation by Death",
                                      c4AD == 0 & c4Dementia == 0 & dies == 0 ~ "Survivor"
                                      ))
act_person_last <- act_visit %>%
  filter(row_number() == max(row_number()))

table(act_person_last$c4ADDementia)

# CASI threshold and dementia diagnosis
act_visit_threshold <- act_visit %>%
  filter(row_number() == max(row_number())) %>%
  mutate(casi_over_threshold = casi_score < 85)
with(act_visit_threshold, table(anydementia, casi_over_threshold, useNA = "always"))
with(act_visit_threshold, table(anyad, casi_over_threshold, useNA = "always"))
# for 251 patient the CASI score is not below threshold but still dementia diagnosis
# those are mostly patients whose casi score lies very close to (but above) the threshold
act_visit_threshold %>%
  filter(row_number() == max(row_number())) %>%
  filter(casi_over_threshold == FALSE & anydementia == 1) %>%
  pull(casi_score) %>%
  summary()

# casi score for survivors and patients who die
act_visit %>%
  ungroup() %>%
  filter(!is.na(casi_score)) %>%
  group_by(!is.na(Age_Death)) %>%
  summarize(mean(casi_score))

act_visit %>%
  filter(row_number() == max(row_number())) %>%
  ungroup() %>%
  group_by(c4dth) %>%
  filter(!is.na(casi_score)) %>%
  summarize(mean(casi_score))

# Covariates --------------------------------------------------------------
#check which variables actually vary over time and which ones don't
nunique <- function(x) length(unique(x))
nunique_summ <- act_visit %>%
  group_by(SUBJECT) %>%
  summarize(across(everything(),
                   .fns=list(nunique=nunique),
                   .names = "{.fn}_{.col}"))



# Age ---------------------------------------------------------------------
# create columns for age at baseline, age at AD / dementia diagnosis / death
act_visit <- act_visit %>%
  group_by(SUBJECT) %>%
  mutate(visitdt_bl = visitdt,
         visitdt_lag = lag(visitdt),
         visitdt_last = max(visitdt),
         age_bl = min(age)) %>%
  mutate(age_new = age_bl + as.numeric(visitdt_bl - min(visitdt_bl)) / 365.25,
         age_death = ifelse(!is.na(Age_Death),
                            age_bl + as.numeric(deathdt - min(visitdt_bl)) / 365.25,
                            NA)) %>%
  mutate(age_ADdem = ifelse(c4Dementia == 1,
                            ((age_bl + as.numeric(visitdt_lag - min(visitdt_bl))/365.25) + age_new) / 2,
                            NA))


# filter out prevalent cases: no prevalent cases
act_visit %>% filter(max(row_number())==1 & anyad ==1)

# range of age (censoring applied)
act_visit$age_new %>% range()
act_visit$age_death %>% na.omit() %>% range()

act_visit %>% filter(c4AD == 1) %>% pull(age_new) %>% range()
act_visit %>% group_by(SUBJECT) %>% filter(row_number() == max(row_number()) & c4AD == 0) %>% pull(age_new) %>% range()


act_visit %>% filter(!is.na(deathdt)) %>% pull(age) %>% range()

# average age for Dementia diagnosis (after imputation)
act_visit %>% filter(c4Dementia == 1) %>% pull(age_new) %>% mean()
# average age for AD diagnosis (after imputation)
act_visit %>% filter(c4AD == 1) %>% pull(age_new) %>% mean()

# Depletion of susceptibles  ----------------------------------------------
act_person_final <- act_visit %>%
  filter(VISIT == max(VISIT))
act_person_final %>% filter(age_death < 85) %>% nrow()
act_person_final %>% filter(age_death >= 85) %>% nrow()

act_person_final %>% filter(age_death < 85 & c4AD == 1) %>% pull(SUBJECT) %>% unique() %>% length()

act_person_final %>% filter(age_death >= 85 & c4AD == 1) %>% pull(SUBJECT) %>% unique() %>% length()


with(act_person_final, prop.table(table(age_death < 85, c4Dementia)))
with(act_person_final, prop.table(table(age_death < 85, c4AD), margin = 1))

# apoe
act_person_final$apoe %>%
  table(useNA = "always") %>%
  prop.table()
act_person_final %>%
  filter(age_death < 85) %>%
  pull(apoe) %>%
  table(useNA = "always") %>%
  prop.table()
act_person_final %>%
  filter(age_death >= 85) %>%
  pull(apoe) %>%
  table(useNA = "always") %>%
  prop.table()

# smoking
act_person_final$smoker_ind %>%
  table(useNA = "always") %>%
  prop.table()
act_person_final %>%
  filter(age_death < 85) %>%
  pull(smoker_ind) %>%
  table(useNA = "always") %>%
  prop.table()
act_person_final %>%
  filter(age_death >= 85) %>%
  pull(smoker_ind) %>%
  table(useNA = "always") %>%
  prop.table()

# Visualize ---------------------------------------------------------------
# split patients into patient groups for better visualization
groups <- act_visit$SUBJECT %>%
  unique() %>%
  split(f = cut(seq(1, N),
                breaks = 115,
                labels = FALSE))

# sample_id <- act_visit$SUBJECT %>%
#   unique() %>%
#   sample(size = 50, replace = FALSE)

for (i in seq(1, 115)) {
  # filter patients in that subgroup
  act_visit_sample <- act_visit %>%
    filter(SUBJECT %in% groups[[i]]) %>%
    filter(casi_is_valid == 1)

  # act_visit_sample$age_new %>% is.na() %>% any()
  # act_visit_sample$casi_score %>% is.na() %>% any()
  # act_visit_sample$c4Dementia %>% is.na() %>% any()
  # act_visit_sample$c4AD %>% is.na() %>% any()
  act_visit_sample$dies

  ggplot(data = act_visit_sample,
         aes(x = age_new,
             y = casi_score,
             group = SUBJECT,
             color = ADDementia_color %>% as.factor())) +
    geom_line() +
    geom_hline(aes(yintercept = 85),
               color = "grey",
               linetype = "dashed") +
    geom_vline(aes(xintercept = 65),
               color = "grey",
               linetype = "dashed") +
    geom_point(size = 3) +
    # geom_point(aes(x = age_ADdem,
    #                y = 65),
    #            shape = 4,
    #            color = "red") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(color = "Dementia / AD / Death?",
         x = "Age in years", y = "CASI score") +
    scale_color_manual(values = c("AD" = "#880808",
                                  "other dementia" = "#6495ED",
                                  "Truncation by Death" = "#636363",
                                  "Survivor"  = "#8A9A5B"))

  ggsave(filename = paste0("ACT/Results/CASI_trajectory_subgroup", i, ".pdf"),
         unit = "cm",
         width = 25, height = 18)
}


# average casi score for AD / dementia patients
AD_patients <- act_visit %>% filter(c4AD == 1) %>% pull(SUBJECT)
Dementia_patients <- act_visit %>% filter(c4Dementia == 1 & c4AD != 1) %>% pull(SUBJECT)

act_AD_av <- act_visit %>%
  filter(SUBJECT %in% AD_patients) %>%
  ungroup() %>%
  group_by(age_new) %>%
  summarize(casi_score_av = mean(na.omit(casi_score))) %>%
  ungroup() %>%
  mutate(type = "AD")

ggplot(data = act_AD_av,
       aes(x = age_new,
           y = casi_score_av)) +
  # geom_line() +
  geom_hline(aes(yintercept = 85),
             color = "grey",
             linetype = "dashed") +
  geom_vline(aes(xintercept = 65),
             color = "grey",
             linetype = "dashed") +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  theme_bw() +
  labs(x = "Age in years",
       y = "CASI score",
       title = "Average CASI score among patients who will be diagnosed with AD")


act_Dementia_av <- act_visit %>%
  filter(SUBJECT %in% Dementia_patients) %>%
  ungroup() %>%
  group_by(age_new) %>%
  summarize(casi_score_av = mean(na.omit(casi_score))) %>%
  ungroup() %>%
  mutate(type = "other dementia")

ggplot(data = act_Dementia_av,
       aes(x = age_new,
           y = casi_score_av)) +
  # geom_line() +
  geom_hline(aes(yintercept = 85),
             color = "grey",
             linetype = "dashed") +
  geom_vline(aes(xintercept = 65),
             color = "grey",
             linetype = "dashed") +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  theme_bw() +
  labs(x = "Age in years",
       y = "CASI score",
       title = "Average CASI score among patients who will be diagnosed with other types of dementia")


# average for healthy
act_visit_av_surv <- act_visit %>%
  filter(!(SUBJECT %in% c(AD_patients, Dementia_patients))) %>%
  ungroup() %>%
  group_by(age_new) %>%
  summarize(casi_score_av = mean(na.omit(casi_score))) %>%
  ungroup() %>%
  mutate(type = "survivors")

ggplot(data = act_visit_av_surv,
       aes(x = age_new,
           y = casi_score_av)) +
  # geom_line() +
  geom_hline(aes(yintercept = 85),
             color = "grey",
             linetype = "dashed") +
  geom_vline(aes(xintercept = 65),
             color = "grey",
             linetype = "dashed") +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  theme_bw() +
  labs(x = "Age in years",
       y = "CASI score",
       title = "Average CASI score among survivors and no AD / no dementia patients")

act_averages <- rbind(act_visit_av_surv,
                  act_AD_av,
                  act_Dementia_av)


ggplot(data = act_averages,
       aes(x = age_new,
           y = casi_score_av,
           group = type,
           color = type)) +
  geom_hline(aes(yintercept = 85),
             color = "grey",
             linetype = "dashed") +
  geom_vline(aes(xintercept = 65),
             color = "grey",
             linetype = "dashed") +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age in years",
       y = "CASI score",
       title = "Average CASI score",
       color = "")

ggsave(filename = "ACT/Results/average_CASI_trajectory.pdf",
       unit = "cm",
       width = 25, height = 18)


# Odds ratio --------------------------------------------------------------
# plot odds ratio for death given AD or no AD over age [in years]
act_odds <- act_visit %>%
  ungroup() %>%
  group_by(age = round(age_new, 0)) %>%
  summarize(p_1_1 = sum(c4dth == 1 & c4AD == 1) / sum(c4AD),
            p_0_1 = sum(c4dth == 0 & c4AD == 1) / sum(c4AD),
            p_1_0 = sum(c4dth & c4AD == 0) / sum(c4AD == 0),
            p_0_0 = sum(c4dth == 0 & c4AD == 0) / sum(c4AD == 0),
            at_risk = n()) %>%
  # mutate(p_1_1 = ifelse(p_1_1 == 0, 0.01, p_1_1),
  #        p_0_1 = ifelse(p_0_1 == 0, 0.01, p_0_1),
  #        p_1_0 = ifelse(p_1_0 == 0, 0.01, p_1_0),
  #        p_0_0 = ifelse(p_0_0 == 0, 0.01, p_0_0)
  #        ) %>%
  mutate(odds = (p_1_1 / p_0_1) / (p_1_0 / p_0_0))

# average odds for > 65 years old
act_odds %>%
  filter(!is.nan(odds) & odds != Inf) %>%
  pull(odds) %>%
  mean()
# average odds for > 85 years old
act_odds %>%
  filter(!is.nan(odds) & odds != Inf & age >= 85) %>%
  pull(odds) %>%
  mean()


# APOEe4 allele -----------------------------------------------------------
# joint odds ratio (like Nevo) for AD and Death by APOE allele
act_joint_odds <- act_visit %>%
  filter(!is.na(apoe)) %>%
  group_by(SUBJECT) %>%
  mutate(c4dth_lag = lag(c4dth, 1),
         c4AD_lag = lag(c4AD)) %>%
  mutate(c4dth_lag = ifelse(is.na(c4dth_lag), 0, c4dth_lag),
         c4AD_lag = ifelse(is.na(c4AD_lag), 0, c4AD_lag)) %>%
  ungroup() %>%
  group_by(age = round(age_new, 0),
           apoe) %>%
  summarize(p_1_1 = sum(c4dth == 1 & c4AD == 1 & c4dth_lag == 0 & c4AD_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            p_0_1 = sum(c4dth == 0 & c4AD == 1 & c4AD_lag == 0 & c4dth_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            p_1_0 = sum(c4dth == 1 & c4AD == 0 & c4AD_lag == 0 & c4dth_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            p_0_0 = sum(c4dth == 0 & c4AD == 0 & c4AD_lag == 0 & c4dth_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            at_risk = n()) %>%
  # mutate(p_1_1 = ifelse(p_1_1 == 0, 0.01, p_1_1),
  #        p_0_1 = ifelse(p_0_1 == 0, 0.01, p_0_1),
  #        p_1_0 = ifelse(p_1_0 == 0, 0.01, p_1_0),
  #        p_0_0 = ifelse(p_0_0 == 0, 0.01, p_0_0)
  #        ) %>%
  mutate(odds = (p_1_1 / p_0_1) / (p_1_0 / p_0_0))

# plot joint odds ratio for AD and death
# for APOE 1 or 0
ggplot(data = act_joint_odds,
       aes(x = age,
           y = odds,
           group = apoe,
           color = apoe %>% as.factor())) +
  geom_point() +
  geom_line(linetype = "dashed") +
  geom_vline(aes(xintercept = 65),
             color = "#4F7942") +
  geom_vline(aes(xintercept = 85),
             color = "#50C878") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age [in full years]",
       title = "Joint odds ratio for death and AD diagnosis",
       y = "Odds ratio",
       color = "APOEe4 allele") +
  scale_color_manual(values = c("0" = "#6082B6",
                                "1" = "#000080"))

ggsave("ACT/Results/joint_odds_ratio.pdf",
       unit = "cm",
       width = 20, height = 15)


# plot odds ratio for death given APOE over age [in years]
act_odds <- act_visit %>%
  filter(!is.na(apoe)) %>%
  ungroup() %>%
  group_by(age = round(age_new, 0)) %>%
  summarize(p_1_1 = sum(c4AD == 1 & apoe == 1) / sum(apoe),
            p_0_1 = sum(c4AD == 0 & apoe == 1) / sum(apoe),
            p_1_0 = sum(c4AD & apoe == 0) / sum(apoe == 0),
            p_0_0 = sum(c4AD == 0 & apoe == 0) / sum(apoe == 0),

            p_1_1_d = sum(c4dth == 1 & apoe == 1) / sum(apoe),
            p_0_1_d = sum(c4dth == 0 & apoe == 1) / sum(apoe),
            p_1_0_d = sum(c4dth & apoe == 0) / sum(apoe == 0),
            p_0_0_d = sum(c4dth == 0 & apoe == 0) / sum(apoe == 0),
            at_risk = n()) %>%
  # mutate(p_1_1 = ifelse(p_1_1 == 0, 0.01, p_1_1),
  #        p_0_1 = ifelse(p_0_1 == 0, 0.01, p_0_1),
  #        p_1_0 = ifelse(p_1_0 == 0, 0.01, p_1_0),
  #        p_0_0 = ifelse(p_0_0 == 0, 0.01, p_0_0)
  #        ) %>%
  mutate(odds_apoe = (p_1_1 / p_0_1) / (p_1_0 / p_0_0),
         odds_death = (p_1_1_d / p_0_1_d) / (p_1_0_d / p_0_0_d))


ggplot(data = act_odds,
       aes(x = age,
           y = odds_apoe,
           color = "AD")) +
  geom_point() +
  geom_line(linetype = "dashed") +
  geom_point(aes(y = odds_death,
                 color = "Death")) +
  geom_line(aes(y = odds_death,
                color = "Death"),
            linetype = "dashed") +
  geom_vline(aes(xintercept = 65),
             linetype = "dashed",
             color = "grey") +
  geom_hline(aes(yintercept = 1),
             linetype = "dashed",
             color = "grey") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age [in full years]",
       y = "Odds ratio",
       title = "Odds ratio for mortality / AD among patients with and without APOE allel",
       color = "Event?")
ggsave("ACT/Results/mortality_APOE_odds_ratio.pdf",
       unit = "cm",
       width = 20, height = 15)


# Smoking -----------------------------------------------------------------
# joint odds ratio (like Nevo) for AD and Death by smoking status
act_joint_odds_smoking <- act_visit %>%
  filter(!is.na(smoker_ind)) %>%
  group_by(SUBJECT) %>%
  mutate(c4dth_lag = lag(c4dth, 1),
         c4AD_lag = lag(c4AD)) %>%
  mutate(c4dth_lag = ifelse(is.na(c4dth_lag), 0, c4dth_lag),
         c4AD_lag = ifelse(is.na(c4AD_lag), 0, c4AD_lag)) %>%
  ungroup() %>%
  group_by(age = round(age_new, 0),
           smoker_ind) %>%
  summarize(p_1_1 = sum(c4dth == 1 & c4AD == 1 & c4dth_lag == 0 & c4AD_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            p_0_1 = sum(c4dth == 0 & c4AD == 1 & c4AD_lag == 0 & c4dth_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            p_1_0 = sum(c4dth == 1 & c4AD == 0 & c4AD_lag == 0 & c4dth_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            p_0_0 = sum(c4dth == 0 & c4AD == 0 & c4AD_lag == 0 & c4dth_lag == 0) / sum(c4AD_lag == 0 & c4dth_lag == 0),
            at_risk = n()) %>%
  # mutate(p_1_1 = ifelse(p_1_1 == 0, 0.01, p_1_1),
  #        p_0_1 = ifelse(p_0_1 == 0, 0.01, p_0_1),
  #        p_1_0 = ifelse(p_1_0 == 0, 0.01, p_1_0),
  #        p_0_0 = ifelse(p_0_0 == 0, 0.01, p_0_0)
  #        ) %>%
  mutate(odds = (p_1_1 / p_0_1) / (p_1_0 / p_0_0))

act_joint_odds_smoking %>%
  filter(age < 85) %>%
  group_by(smoker_ind) %>%
  summarize(mean(odds))

# plot joint odds ratio for AD and death
# for APOE 1 or 0
ggplot(data = act_joint_odds_smoking,
       aes(x = age,
           y = odds,
           group = smoker_ind,
           color = smoker_ind %>% as.factor())) +
  geom_point() +
  geom_line(linetype = "dashed") +
  geom_vline(aes(xintercept = 65),
             color = "#4F7942") +
  geom_vline(aes(xintercept = 85),
             color = "#50C878") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age [in full years]",
       title = "Joint odds ratio for death and AD diagnosis",
       y = "Odds ratio",
       color = "Smoking status",
       caption = "never smoker N = 2814 \n former / current smoker N = 2934") +
  scale_color_manual(values = c("never smoker" = "#6082B6",
                                "former / current smoker" = "#000080"))
ggsave("ACT/Results/mortality_smoking_odds_ratio.pdf",
       unit = "cm",
       width = 20, height = 15)



