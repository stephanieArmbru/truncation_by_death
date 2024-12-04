
#two goals

#1. make a simple dataset with one row per subject with
  # indicator for ad/dementia onset
  # interval for ad/dementia onset
  # midpoint-imputed time of ad/dementia onset
  # indicator for death
  # time of death
  # in the absence of death, administrative censoring time
    # not clear in ACT whether death monitored for all participants, so may just use time of last visit


library(haven)
library(dplyr)
library(tidyr)
library(SemiCompRisksFreq)

#Load ACT data
ACTpath <- "/Users/htr2/Documents/ACT/"
# ACTscripts <- "/Users/reederh/Dropbox/Harrison Reeder/ACT/Scripts/"
ACTinput <- paste0(ACTpath,"Input/")
ACTtemppath <- paste0(ACTpath,"Temp/")
ACToutpath <- paste0(ACTpath,"Output/")


#read in full visit-level data
act_visit <- haven::read_sas(data_file = paste0(ACTinput,"bienvisit.sas7bdat"))

#make table just of baseline visits to have person-level info
act_person <- act_visit %>% group_by(SUBJECT) %>% filter(VISIT==0)

#look at all the variables
summary(act_person)
summary(act_visit)

#it appears the data were pulled on or around 2020-03-05, so we can use that as the
#'administrative censoring' date
max(act_visit$visitdt)
max(act_visit$deathdt,na.rm=TRUE)

#look at ad/dementia diagnoses by visit
with(act_visit, table(VISIT, anydementia|anyad, useNA = "always"))
#look at correspondence of dementia flag and dsmivdx diagnostic categorical
with(act_visit, table(dsmivdx, anydementia, useNA = "always"))
#make sure everyone with a dementia onset date falls into a dementia category
with(act_visit, table(is.na(onsetdate), dsmivdx, useNA = "always"))
#make sure everyone with a dementia age falls into a dementia category
with(act_visit, table(is.na(Age_Onset), dsmivdx, useNA = "always"))


#look at correspondence of ad flag and nindx diagnostic categorical
with(act_visit, table(nindx, anyad, useNA = "always"))
#look at correspondence of dementia and nindx diagnostic categoricals (not complete overlap!)
with(act_visit, table(nindx, dsmivdx, useNA = "always"))

#check to make sure there were no visits that happened after death (there is one!!)
# View(act_visit[which(act_visit$visitdt > act_visit$deathdt),])
# View(act_visit %>% filter(SUBJECT == 154130))

#check which variables actually vary over time and which ones don't
nunique <- function(x) length(unique(x))
nunique_summ <- act_visit %>% group_by(SUBJECT) %>% summarize(across(everything(),
                                                     .fns=list(nunique=nunique), .names = "{.fn}_{.col}"))
#if this is greater than 1, then at least one subject has multiple unique values of the variable at different visits
#and it is therefore time varying in some way...
apply(nunique_summ,MARGIN = 2,FUN = max)

# act_visit %>% group_by(SUBJECT) %>% 
#   mutate(nunique = nunique(angina)) %>% 
#   filter(nunique > 1) %>% 
#   select(SUBJECT,VISIT,angina,Angina_Yr) %>% View
# act_visit %>% select(SUBJECT,age) %>% group_by(SUBJECT,age) %>% 
#   mutate(count_unique_age = n()) %>% filter(count_unique_age>1)


#make dataset just saying which visits have first AD/dementia onset
act_addemage <- act_visit %>% 
  #for now, consider first observation with observed AD or dementia by any criterion
  filter(anydementia == 1 | anyad == 1) %>% 
  group_by(SUBJECT) %>% filter(row_number()==1) %>% mutate(first_addem=1) %>% 
  select(SUBJECT,VISIT,first_addem) %>% ungroup

#add on number and date of previous visit (to get 'interval' of diagnosis)
act_addemage2 <- act_visit %>% left_join(act_addemage) %>%
  group_by(SUBJECT) %>%
  mutate(
    visit_first_addem=VISIT,
    visit_first_addem_prev=lag(VISIT),
    visitdt_first_addem=visitdt,
    visitdt_first_addem_prev=lag(visitdt)) %>%
  select(SUBJECT,visit_first_addem,visit_first_addem_prev,
         visitdt_first_addem,visitdt_first_addem_prev,
         first_addem) %>% filter(first_addem==1)

#make dataset just with date and age at last visit
act_lastage <- act_visit %>% group_by(SUBJECT) %>% 
  filter(row_number()==max(row_number())) %>% ungroup %>%
  mutate(visit_last=VISIT,visitdt_last=visitdt) %>% 
  select(SUBJECT,visit_last,visitdt_last)

#now, pull together a cleaned person-level dataset with these additional time variables
act_person_cleaned <- act_person %>% rename(visitdt_bl=visitdt, age_bl = age, age_death_old = Age_Death) %>%
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
  select(-gender,-VISIT) %>%
  #look at all the time-varying variables and relabel them with "_bl" suffix
  #to clarify that they're baseline variables
  rename(martial_bl = marital, bmi_bl=bmi, bmi4_bl=bmi4,
         smoke_bl=smoke,exercise_reg_bl=exercise_reg,
         Gait_able_bl=Gait_able,
         Gait_Aid_bl=Gait_Aid,Gait_Time_bl=Gait_Time,
         Chair_Able_bl=Chair_Able,Chair_Time_bl=Chair_Time,
         grip_able_bl=grip_able,grip_strength_bl=grip_strength,
         CESD_Score_bl=CESD_Score,CESD_Flag_bl=CESD_Flag,
         adl_sum_bl=adl_sum,adl_flag_bl=adl_flag,iadl_sum_bl=iadl_sum,iadl_flag_bl=iadl_flag,
         COPD_bl=COPD,Asthma_bl=Asthma,pneumonia_bl=pneumonia,COPD_Yr_bl=COPD_Yr,Asthma_Yr_bl=Asthma_Yr,
         Hypertension_bl=Hypertension,Hypertension_Yr_bl=Hypertension_Yr,CHF_bl=CHF,CHF_Yr_bl=CHF_Yr,
         Stroke_bl=Stroke,TIA_bl=TIA,CEA_bl=CEA,Stroke_Yr_bl=Stroke_Yr,TIA_Yr_bl=TIA_Yr,CEA_Yr_bl=CEA_Yr,
         MI_bl=MI,angina_bl=angina,CABG_bl=CABG,angio_bl=angio,
         MI_Yr_bl=MI_Yr,Angina_Yr_bl=Angina_Yr,CABG_Yr_bl=CABG_Yr,Angio_Yr_bl=Angio_Yr,
         casi_valid_bl=casi_valid,CASI_SC_bl=CASI_SC, CASISHRT_bl=CASISHRT,
         dsmivdx_bl=dsmivdx,nindx_bl=nindx,
         anyad_bl=anyad,anydementia_bl=anydementia,onsetdate_bl=onsetdate,Age_Onset_bl=Age_Onset
         ) %>%
  mutate(CASI_SC_corrected_bl = ifelse(CASI_SC_bl > 100, NA, CASI_SC_bl))

#generate final semi-competing risk outcomes
#yL: age at baseline (setting age 65=0)
#y1: age at ad/dementia onset (setting age 65=0)
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



summary((act_person_cleaned$y1 - act_person_cleaned$yL))

summary((act_person_cleaned$y2 - act_person_cleaned$y1)[act_person_cleaned$delta1==1])
summary((act_person_cleaned$y2 - act_person_cleaned$y1)[act_person_cleaned$delta1==0])
summary((act_person_cleaned$y2 - act_person_cleaned$y1)[act_person_cleaned$delta1==0])

fit_unadj <- FreqID_HReg2(Formula(yL | y1 + delta1 | y2 + delta2 ~ 1 | 1 | 1), 
                #NOTE THAT SOME PEOPLE HAVE ONE VISIT (AGE_BL==AGE_LAST) AND THAT THROWS THINGS OFF, SO THEY NEED TO BE REMOVED
                data = act_person_cleaned %>% filter(age_last != age_bl), 
                hazard = "bs", frailty = TRUE, model = "semi-markov")
plot(predict(fit_unadj),plot.est = "Surv")

