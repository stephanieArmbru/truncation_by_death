################################### ACT data ###################################
############################### Example Patients ###############################

##### Stephanie Armbruster
##### September 2024



# Library -----------------------------------------------------------------
library(tidyverse)
library(ggpubr)



# Shared frailty ID model -------------------------------------------------
# Semicompeting risks -----------------------------------------------------
#make dataset just saying which visits have first AD/dementia onset
act_addemage <- act_visit %>%
  #for now, consider first observation with observed AD or dementia by any criterion
  filter(anydementia == 1 | anyad == 1) %>%
  group_by(SUBJECT) %>%
  filter(row_number() == 1) %>%
  mutate(first_addem = 1) %>%
  dplyr::select(SUBJECT, VISIT, first_addem) %>%
  ungroup()

#add on number and date of previous visit (to get 'interval' of diagnosis)
act_addemage2 <- act_visit %>%
  left_join(act_addemage) %>%
  group_by(SUBJECT) %>%
  mutate(
    visit_first_addem=VISIT,
    visit_first_addem_prev=lag(VISIT),
    visitdt_first_addem=visitdt,
    visitdt_first_addem_prev=lag(visitdt)) %>%
  dplyr::select(SUBJECT,visit_first_addem,visit_first_addem_prev,
                visitdt_first_addem,visitdt_first_addem_prev,
                first_addem) %>%
  filter(first_addem==1)

#make dataset just with date and age at last visit
act_lastage <- act_visit %>%
  group_by(SUBJECT) %>%
  filter(row_number()==max(row_number())) %>%
  ungroup() %>%
  mutate(visit_last=VISIT,visitdt_last=visitdt) %>%
  dplyr::select(SUBJECT,visit_last,visitdt_last)

#now, pull together a cleaned person-level dataset with these additional time variables
act_person_cleaned <- act_person %>%
  rename(visitdt_bl=visitdt,
         age_bl = age,
         age_death_old = Age_Death) %>%
  left_join(act_addemage2) %>%
  left_join(act_lastage) %>%
  #make slightly more 'precise' age variables based on baseline age and then visit dates
  #because currently they are rounded to the nearest integer and that can cause some trouble
  #obviously it's a bit imprecise because the age_bl is itself rounded to the nearest integer,
  #but at least it's more precise relative to that...
  mutate(
    age_first_addem = age_bl + as.numeric(visitdt_first_addem - visitdt_bl)/365.25,
    age_first_addem_prev = age_bl + as.numeric(visitdt_first_addem_prev - visitdt_bl)/365.25,
    age_first_addem_mid = (age_first_addem_prev + age_first_addem) / 2,
    age_death = age_bl + as.numeric(deathdt - visitdt_bl)/365.25,
    pulldt = as.Date("2020-03-05"),
    age_last = age_bl + as.numeric(visitdt_last - visitdt_bl)/365.25,
    age_pulled = age_bl + as.numeric(pulldt - visitdt_bl)/365.25) %>%
  #drop one weird study visit which appears to have occurred after death...
  filter(!(SUBJECT == 154130 & VISIT==6)) %>%
  mutate(apoe_raw = ifelse(apoe_raw=="",NA,apoe_raw)) %>%
  separate(col=apoe_raw,into=c("apoe1","apoe2")) %>%
  mutate(female = gender - 1, #create more interpretable variable name
         apoe4_any = as.numeric(pmax(apoe1,apoe2)==4), #create meaningful APOE variables
         apoe4_ct = as.numeric(apoe1==4) + as.numeric(apoe2==4)) %>%
  dplyr::select(-gender,-VISIT) %>%
  #look at all the time-varying variables and relabel them with "_bl" suffix
  #to clarify that they're baseline variables
  rename(martial_bl = marital,
         bmi_bl=bmi,
         bmi4_bl=bmi4,
         smoke_bl=smoke,
         exercise_reg_bl=exercise_reg,

         Gait_able_bl=Gait_able,
         Gait_Aid_bl=Gait_Aid,
         Gait_Time_bl=Gait_Time,

         Chair_Able_bl=Chair_Able,
         Chair_Time_bl=Chair_Time,

         grip_able_bl=grip_able,
         grip_strength_bl=grip_strength,

         CESD_Score_bl=CESD_Score,
         CESD_Flag_bl=CESD_Flag,

         adl_sum_bl=adl_sum,
         adl_flag_bl=adl_flag,
         iadl_sum_bl=iadl_sum,
         iadl_flag_bl=iadl_flag,

         COPD_bl=COPD,
         Asthma_bl=Asthma,
         pneumonia_bl=pneumonia,
         COPD_Yr_bl=COPD_Yr,
         Asthma_Yr_bl=Asthma_Yr,
         Hypertension_bl=Hypertension,
         Hypertension_Yr_bl=Hypertension_Yr,
         CHF_bl=CHF,
         CHF_Yr_bl=CHF_Yr,
         Stroke_bl=Stroke,
         TIA_bl=TIA,
         CEA_bl=CEA,
         Stroke_Yr_bl=Stroke_Yr,
         TIA_Yr_bl=TIA_Yr,
         CEA_Yr_bl=CEA_Yr,
         MI_bl=MI,
         angina_bl=angina,
         CABG_bl=CABG,
         angio_bl=angio,
         MI_Yr_bl=MI_Yr,
         Angina_Yr_bl=Angina_Yr,
         CABG_Yr_bl=CABG_Yr,
         Angio_Yr_bl=Angio_Yr,
         casi_valid_bl=casi_valid,
         CASI_SC_bl=CASI_SC,
         CASISHRT_bl=CASISHRT,
         dsmivdx_bl=dsmivdx,
         nindx_bl=nindx,
         anyad_bl=anyad,
         anydementia_bl=anydementia,
         onsetdate_bl=onsetdate,
         Age_Onset_bl=Age_Onset
  ) %>%
  mutate(smoker_bl = case_when(smoke_bl == 0 ~ 0,
                               smoke_bl == 1 | smoke_bl == 2 ~ 1,
                               .default = NA),
         CASI_SC_corrected_bl = ifelse(CASI_SC_bl > 100,
                                       NA,
                                       CASI_SC_bl))

#generate final semi-competing risk outcomes
#yL: age at baseline (setting age 65=0)
#y1: age at ad/dementia onset (setting age 65=0) (middle of interval in which AD/dem must developed as it is intially detected)
#y2: age at death or last study visit (setting age 65=0)
#delta1: indicator for observed ad/dementia onset
#delta2: indicator for observed death

act_person_cleaned$yL <- act_person_cleaned$age_bl - 65
act_person_cleaned$y1 <- case_when(
  !is.na(act_person_cleaned$age_first_addem_mid) ~ act_person_cleaned$age_first_addem_mid,
  !is.na(act_person_cleaned$age_death) ~ act_person_cleaned$age_death,
  #for now, censor based on last visit rather than administrative date, ...
  TRUE ~ act_person_cleaned$age_last) - 65

act_person_cleaned$delta1 <- as.numeric(!is.na(act_person_cleaned$age_first_addem))

act_person_cleaned$y2 <- case_when(
  !is.na(act_person_cleaned$age_death) ~ act_person_cleaned$age_death,
  #for now, censor based on last visit rather than administrative date, ...
  TRUE ~ act_person_cleaned$age_last) - 65
act_person_cleaned$delta2 <- !is.na(act_person_cleaned$age_death)

# check transformed covariates
summary((act_person_cleaned$y1 - act_person_cleaned$yL))

summary((act_person_cleaned$y2 - act_person_cleaned$y1)[act_person_cleaned$delta1==1])
summary((act_person_cleaned$y2 - act_person_cleaned$y1)[act_person_cleaned$delta1==0])
summary((act_person_cleaned$y2 - act_person_cleaned$y1)[act_person_cleaned$delta1==0])

# baseline covariates
act_person_cleaned$smoker_bl %>% table(useNA = "always")
act_person_cleaned$smoke_bl %>% table(useNA = "always")
# N=14 patients have missing smoking data


# no missing data for gender
act_person_cleaned$female %>% table(useNA = "always")

# N = 13 patients have missing exercise data
act_person_cleaned$exercise_reg_bl %>% table(useNA = "always")

# N = 9 patients have missing hispanic data
act_person_cleaned$hispanic %>% table(useNA = "always")

# N = 1077 patients have no genetic data
act_person_cleaned$apoe %>% table(useNA = "always")

# N = 3413 and N = 3616 patients missing
act_person_cleaned$FamHx_Dx %>% table(useNA = "always")
act_person_cleaned$FamHx_AD %>% table(useNA = "always")

# martial level for all patients observed
act_person_cleaned$martial_bl %>% table(useNA = "always")


# fit frailty illness death model
# B-spline hazard function
# Semi-Markov
fit_unadj <- FreqID_HReg2(Formula(yL | y1 + delta1 | y2 + delta2 ~ 1 | 1 | 1),
                          #NOTE THAT SOME PEOPLE HAVE ONE VISIT (AGE_BL==AGE_LAST) AND THAT THROWS THINGS OFF, SO THEY NEED TO BE REMOVED
                          data = act_person_cleaned %>% filter(age_last != age_bl),
                          hazard = "bs",
                          frailty = TRUE,
                          model = "semi-markov")

fit_unadj %>% summary()
# plot survival functions for shared frailty Illness-Death model
plot(predict(fit_unadj),
     plot.est = "Surv")

# plot hazard for shared frailty Illness-Death model
plot(predict(fit_unadj),
     plot.est="Haz")

# complete case analysis
act_person_cleaned_subset <- act_person_cleaned %>%
  filter(!is.na(smoker_bl) & !is.na(female) & !is.na(exercise_reg_bl))

# excluding N=26 patients
N - act_person_cleaned_subset %>% nrow()

# with covariates
fit_cov <- FreqID_HReg2(Formula(yL | y1 + delta1 | y2 + delta2 ~ smoker_bl + female + exercise_reg_bl | smoker_bl + female + exercise_reg_bl | smoker_bl + female + exercise_reg_bl),
                        #NOTE THAT SOME PEOPLE HAVE ONE VISIT (AGE_BL==AGE_LAST) AND THAT THROWS THINGS OFF, SO THEY NEED TO BE REMOVED
                        data = act_person_cleaned_subset %>%
                          filter(age_last != age_bl),
                        hazard = "bs",
                        frailty = TRUE,
                        model = "semi-markov")

fit_cov %>% summary()



# Example risk profiles ---------------------------------------------------
# example patients
# patient 1: male smoker, who does not exercise
# patient 2: female non-smoker, who exercises regularly

example_bl <- data.frame(smoker_bl = c(1, 0),
                         female = c(0, 1),
                         exercise_reg_bl = c(0, 1)) %>%
  as.matrix()

# predict survival probability for patients 1 and 2
plot(predict(fit_cov,
             x1new = example_bl,
             x2new = example_bl,
             x3new = example_bl),
     plot.est = "Surv")

# time frame to predict
tseq_pred <- seq(0, 35)


# predict risk profiles
# baseline function for conditional hazard
h3_b_fun <- function(x) splines2::bSpline(x,
                                          knots = c(0, 20),
                                          Boundary.knots = c(0, Inf),
                                          degree = 0,
                                          intercept = FALSE)

temp_pred <-
  SemiCompRisksFreq:::pred_risk_ID(tseq = tseq_pred,
                                   para = fit_cov$estimate,
                                   # baseline input
                                   x1new = example_bl,
                                   x2new = example_bl,
                                   x3new = example_bl,

                                   frailty = fit_cov$frailty,
                                   model = fit_cov$model,

                                   nP0 = fit_cov$nP0,
                                   nP = fit_cov$nP,
                                   p3tv = 0,
                                   h3tv_basis_func = h3_b_fun, #this is omitted if no t1 effect in h3, we just need to specify the nature of the effect and how many parameters correspond to it
                                   hazard = fit_cov$hazard,
                                   knots_list = fit_cov$knots_list,
                                   n_quad = fit_cov$n_quad,
                                   quad_method = fit_cov$quad_method,
                                   Finv = fit_cov$Finv,
                                   alpha = 0.05)


length_tseq <- length(tseq_pred)




# Visualization -----------------------------------------------------------
plot_frame_1 <- data.frame(
  Time=rep(tseq_pred,4),
  Probability = c(temp_pred$p_neither[,1],
                  temp_pred$p_term_only[,1],
                  temp_pred$p_both[,1],
                  temp_pred$p_nonterm_only[,1]),
  Outcome=factor(x = c(rep("Neither",length_tseq),
                       rep("Death",length_tseq),
                       rep("AD + Death",length_tseq),
                       rep("AD",length_tseq)),
                 levels=c("Neither","Death","AD + Death","AD"))
)

plot_frame_2 <- data.frame(
  Time=rep(tseq_pred,4),
  Probability = c(temp_pred$p_neither[,2],
                  temp_pred$p_term_only[,2],
                  temp_pred$p_both[,2],
                  temp_pred$p_nonterm_only[,2]),
  Outcome=factor(x = c(rep("Neither",length_tseq),
                       rep("Death",length_tseq),
                       rep("AD + Death",length_tseq),
                       rep("AD",length_tseq)),
                 levels=c("Neither","Death","AD + Death","AD"))
)

# risk profile for patient 1
rp_1 <- ggplot(plot_frame_1,
       aes(x=65 + Time, y=Probability)) +
  geom_area(aes(colour=Outcome,
                fill=Outcome)) +
  scale_color_manual(values = four_color_cb)+
  scale_fill_manual(values = four_color_cb) +
  theme_bw() +
  labs(title = "Patient 1: male smoker, never exercise",
       x = "Age [years]")

# risk profile for patient 2
rp_2 <- ggplot(plot_frame_2,
       aes(x=65 + Time, y=Probability)) +
  geom_area(aes(colour=Outcome,
                fill=Outcome)) +
  scale_color_manual(values = four_color_cb)+
  scale_fill_manual(values = four_color_cb) +
  theme_bw() +
  labs(title = "Patient 2: female never smoker, always exercise",
       x = "Age [years]")



ggarrange(rp_1, rp_2,
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = "bottom")


ggsave(filename = "ACT/Results/risk_profiles_examples.pdf",
       width = 25, height = 13, unit = "cm")












