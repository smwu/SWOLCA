#===========================================================
## Data Preparation and Cleaning
## Programmer: SM Wu   
## Data: NHANES 2015-2018  
##   Adult, low-income, women, not breastfeeding or pregnant
## Date Updated: 2023/07/12
#===========================================================

# Supporting Libraries
library(haven)
library(tidyverse)
library(survey)
library(psych)
library(readxl)
library(gmodels)
library(ggpubr)

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Github/wsOFMM/"
data_dir <- "Data/"

#============= Load in dietary and demographic data ============================
food_terts <- read.csv(paste0(wd, data_dir, "nhanes1518_tert_lowF.csv"))

# old_data <- read.csv(paste0(wd, data_dir, "nhanesallage_frscores1118.csv"))

#============= Participants on food stamps or WIC ==============================

# Read in nhanes food security questionnaires
# FSQ165: hh food stamp benefit ever received: 1=y, 2=n, 7=refused, 9=don't know
# FSQ012 (nested): hh received SNAP/food stamp benefits in the last 12 months
# FSQ162: hh received WIC benefit in the last 12 months. Restricted to those 
# eligible for WIC
FED_15 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/FSQ_I.XPT', 
                  col_select = c(SEQN, FSQ012, FSQ162, FSQ165))
FED_17 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/FSQ_J.XPT', 
                  col_select = c(SEQN, FSQ012, FSQ162, FSQ165))
fed_all <- rbind(FED_15, FED_17) 

# Account for nesting: If no fs ever received, then none in last 12 months
fed_all$FSQ012[fed_all$FSQ165 == 2] <- 0
# Create FED_ASSIST variable, equal to 1 if participant on food stamps or WIC
fed_all$FED_ASSIST <- ifelse(fed_all$FSQ012 == 1 | fed_all$FSQ162 == 1, 1, 0)
# Merge in fed assist data
food_fed <- merge(x = food_terts, y = fed_all, by = "SEQN", 
                  all.x = TRUE, all.y = FALSE)
# check missingness
temp <- food_fed[is.na(food_fed$FED_ASSIST), 53:56]  # 351 obs
# Turn missing into no's
food_fed$FED_ASSIST[is.na(food_fed$FED_ASSIST)] <- 0

#==================== Smoking ==================================================

## Smoking Questionnaires (never, former, current, second-hand)
# SMQ020: Smoked at least 100 cigarettes in life
# SMQ040 (nested): Now smoke cigarettes
# SMQ621 (nested): Cigarettes smoked in entire life
# SMQ050Q (nested): How long since quit smoking cigarettes
smo_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SMQ_I.XPT', 
                   col_select = c(SEQN, SMQ020, SMQ621, SMQ050Q, SMQ040))
smo_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SMQ_J.XPT', 
                   col_select = c(SEQN, SMQ020, SMQ621, SMQ050Q, SMQ040))
smo_all <- rbind(smo_15, smo_17)
# According to CDC Tobacco Glossary, define current smoker as having smoked 100 
# cigarettes in lifetime and currently smokes cigarettes
smo_all$smoker <- NA
smo_all$smoker[smo_all$SMQ040 %in% c(1, 2)] <- 1  # current smoker
smo_all$smoker[smo_all$SMQ020 == 2] <- 0  # did not smoke 100 cig in life
smo_all$smoker[smo_all$SMQ020 == 1 & smo_all$SMQ040 == 3] <- 0  # former smoker

# Merge in smoking data
food_fed_smo <- merge(x = food_fed, y = smo_all, by = "SEQN", 
                  all.x = TRUE, all.y = FALSE)
# check missingness
temp <- food_fed_smo[is.na(food_fed_smo$smoker), -c(1:56)]  # 2 obs
# Turn missing into no's
food_fed_smo$smoker[is.na(food_fed_smo$smoker)] <- 0

#==================== Physical activity ========================================

###Physical activity questionnaires: no missingness or refusals!
# PAQ610 (nested): Number of days vigorous work per week
# PAD615 (nested): Minutes vigorous-intensity work per day
# PAQ625 (nested): Number of days moderate work per week
# PAD630 (nested): Minutes moderate-intensity work per day
# PAQ640 (nested): Number of days walk or bicycle per week
# PAD645 (nested): Minutes walk/bicycle for transportation per day
# PAQ655 (nested): Days vigorous recreational activities per week
# PAD660 (nested): Minutes vigorous recreational activities per day
# PAQ670 (nested): Days moderate recreational activities per week
# PAD675 (nested): Minutes moderate recreational activities per day
paq_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/PAQ_I.XPT',
                   col_select = c(SEQN, PAQ610, PAD615, PAQ625, PAD630, PAQ640, 
                                  PAD645, PAQ655, PAD660, PAQ670, PAD675))
paq_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/PAQ_J.XPT', 
                   col_select = c(SEQN, PAQ610, PAD615, PAQ625, PAD630, PAQ640, 
                                  PAD645, PAQ655, PAD660, PAQ670, PAD675))
paq_all <- rbind(paq_15, paq_17)

# Physically active = 1 if >=150 mins of moderate or vigorous exercise per week
paq_all$Mins_Active <- apply(paq_all[,-1], 1, function(x) 
  sum(c(x[1]*x[2], x[3]*x[4], x[5]*x[6], x[7]*x[8], x[9]*x[10]), na.rm = TRUE))
paq_all$Phys_Active <- ifelse(paq_all$Mins_Active >= 150, 1, 0)

# Merge in physical activity data
food_fed_smo_paq <- merge(x = food_fed_smo, 
                          y = paq_all[, c("SEQN", "Mins_Active", "Phys_Active")], 
                          by = "SEQN", all.x = TRUE, all.y = FALSE)
# check missingness
temp <- food_fed_smo_paq[is.na(food_fed_smo_paq$Phys_Active), -c(1:60)]  # 0 obs


#======================= Hypertension ==========================================

## Blood pressure readings from mobile examination center (MEC) data. 
# Systolic (BPXSY) and diastolic (BPXDI) three consecutive readings, with fourth
# if necesssary
bpx_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BPX_I.XPT', 
                     col_select = c(SEQN, BPXSY1, BPXDI1, BPXSY2, BPXDI2, 
                                    BPXSY3, BPXDI3, BPXSY4, BPXDI4))
bpx_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPX_J.XPT', 
                     col_select = c(SEQN, BPXSY1, BPXDI1, BPXSY2, BPXDI2,
                                    BPXSY3, BPXDI3, BPXSY4, BPXDI4))
bpx_all <- rbind(bpx_15, bpx_17)
# Take average over four readings, then assign hypertension to those above 130/80
# average systolic bp
bpx_all$SBP_avg <- rowMeans(cbind(bpx_all$BPXSY1, bpx_all$BPXSY2, bpx_all$BPXSY3,
                                  bpx_all$BPXSY4), na.rm=TRUE)
# average diastolic bp
bpx_all$DBP_avg <- rowMeans(cbind(bpx_all$BPXDI1, bpx_all$BPXDI2, bpx_all$BPXDI3,
                                 bpx_all$BPXDI4), na.rm=TRUE)
# flag hypertension per AHA guidelines: https://www.heart.org/en/health-topics/high-blood-pressure
bpx_all$BPX_flag <- ifelse(bpx_all$SBP_avg > 130 | bpx_all$DBP_avg > 80, 1, 0)  

## Blood pressure questionnaire: none refused or missing! 10 don't know
# BPQ020: Ever told you had high blood pressure
# BPQ030: Told had high blood presuure 2+ times
# BPQ050A: Now taking prescribed medicine for HBP
# Assign hypertension to those currently taking blood pressure medication
# or those told they had hypertension 2+ times
bpq_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BPQ_I.XPT', 
                   col_select = c(SEQN, BPQ020, BPQ030, BPQ050A))
bpq_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPQ_J.XPT', 
                   col_select = c(SEQN, BPQ020, BPQ030, BPQ050A))
bpq_all <- rbind(bpq_15, bpq_17)
bpq_all$BPQ_flag <- ifelse(bpq_all$BPQ030 == 1, 1, 0)  # told HPB 2+ times
bpq_all$BPQ_flag[bpq_all$BPQ050A == 1] <- 1  # taking medication
bpq_all$BPQ_flag[bpq_all$BPQ020 == 2] <- 0

# Merge in hypertension data
food_fed_smo_paq_bp <- food_fed_smo_paq %>% 
  left_join(bpx_all %>% dplyr::select(SEQN, SBP_avg, DBP_avg, BPX_flag), by = "SEQN") %>%
  left_join(bpq_all, by = "SEQN")
# Define hypertension flag as flagged from readings or from questionnaire
food_fed_smo_paq_bp$BP_flag <- ifelse(food_fed_smo_paq_bp$BPX_flag == 1 | 
                                      food_fed_smo_paq_bp$BPQ_flag == 1, 1, 0)
# check missingness
temp <- food_fed_smo_paq_bp[is.na(food_fed_smo_paq_bp$BP_flag), -c(1:62)]  # 35 obs
# Turn missing into no's
food_fed_smo_paq_bp$BP_flag[is.na(food_fed_smo_paq_bp$BP_flag)] <- 0


#===================== Save dataset ============================================

write.csv(food_fed_smo_paq_bp,
          paste0(wd, data_dir, "nhanes1518_adult_low_f_12jul2023.csv"))




#===================== Miscellaneous old code ==================================

# #============= Demographic variables ===========================================
# # gender, age, race, education, weight, masked variance pseudo, 
# # income poverty ratio, federal assistance
# DEMO_15 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT', 
#                    col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,
#                                   INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
# DEMO_17 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT', 
#                    col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,
#                                   INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
# # Restrict to adults age 20 or older
# DEMO_15 <- DEMO_15[ which(DEMO_15$RIDAGEYR >= 20), ]
# DEMO_17 <- DEMO_17[ which(DEMO_17$RIDAGEYR >= 20), ]
# demo_temp <- rbind(DEMO_15, DEMO_17) 
# demo_temp$RIAGENDR <- factor(demo_temp$RIAGENDR, levels=c(1, 2), 
#                              labels=c("Male", "Female"))
# demo_temp$RIDRETH3 <- factor(demo_temp$RIDRETH3, levels=c(1:4, 6, 7), 
#                              labels = c("Mexican","Other Hispanic", "NH White", 
#                                         "NH Black", "NH Asian", "Other/Mixed"))
# demo_temp$EDU <- ifelse(demo_temp$DMDEDUC2 < 3, "less than HS",
#                         ifelse(demo_temp$DMDEDUC2 == 3, "HS/GED", 
#                                ifelse(demo_temp$DMDEDUC2 > 3,
#                                       "At least some college", "NA")))
# demo_temp$WTINT4YR <- demo_temp$WTINT2YR/2
# demo_temp$WTMEC4YR <- demo_temp$WTMEC2YR/2
# 

### SMOKING
# smo_fin <- 
#   mutate(smo_all, Smoking_Status = case_when(
#     SMQ621 == 1 ~ "Never Smoked",
#     SMQ040 == 3 ~ "Former Smoker",
#     SMQ020 == 1 | SMQ621 %in% c(2:8) ~ "Current Smoker",
#     SMQ020 == 2 ~ "Light Smoker"
#   )) 
# View(smo_fin)
# smo_fin$smoke_num <- ifelse(smo_fin$Smoking_Status=="Never Smoked" | 
#                               smo_fin$Smoking_Status=="Light Smoker", 1, 
#                             ifelse(smo_fin$Smoking_Status=="Former Smoker", 2, 3))
# smo_fin$smoke_num <- factor(smo_fin$smoke_num, c(1:3),
#                             labels=c("Not Regular", "Former", "Current"))
# smo_keep <- smo_fin[, c("SEQN", "smoke_num")]
# 
# ## Household Smoking Data
# hh_smo_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SMQFAM_I.XPT', 
#                       col_select = c(SEQN, SMD470))
# hh_smo_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SMQFAM_J.XPT', 
#                       col_select = c(SEQN, SMD470))
# hh_smo_all <- rbind(hh_smo_15, hh_smo_17)
# hh_smo_all$Household_Smoker <- ifelse(hh_smo_all$SMD470 == 0, 1, 2)
# hh_smo_all$Household_Smoker <- factor(hh_smo_all$Household_Smoker, c(1,2), 
#                                       labels=c("No","Yes"))
# hh_smo <- hh_smo_all[,c("SEQN","Household_Smoker")]
# # Add household smoking to demo dataset
# demo_smo_hh <- merge(x = demo_smo, y = hh_smo, by = "SEQN", 
#                      all.x = TRUE, all.y = FALSE)

### Physical Activity
# #NA + NA = NA , NA + n = n
# sum1 <- function(x){
#   if(all(is.na(x))){
#     return(NA)
#   } else{
#     return(sum(x, na.rm = TRUE))
#   }
# }
# paq_all <- paq_all %>%
#   mutate(Min_Active1 = case_when (!is.na(PAQ610) & !is.na(PAD615) ~ PAQ610 * PAD615))
# paq_all <- paq_all %>%
#   mutate(Min_Active2 = case_when (!is.na(PAQ640) & !is.na(PAD645) ~ PAQ640 * PAD645))
# paq_all <- paq_all %>%
#   mutate(Min_Active3 = case_when (!is.na(PAQ655) & !is.na(PAD660) ~ PAQ655 * PAD660))
# paq_all$Min_Active <- apply(X = paq_all[, 8:10], MARGIN = 1, FUN = sum1)
# paq_all <- paq_all %>%
#   mutate(Phys_Active = case_when ( 
#     Min_Active >= 150 ~ "Active",
#     Min_Active < 150 ~ "Not Active"))
# paq_all <- paq_all[, c("SEQN", "Phys_Active")]

### Hypertension
# bpq_all$bpq_num <- ifelse(bpq_all$BPQ030 == 1, 1, 2)
# # Nested questionnaire so those with missing converted to No's
# bpq_all$bpq_num[is.na(bpq_all$bpq_num)] <- 2  
# bpq_all$bpq_num <- factor(bpq_all$bpq_num, levels = c(1, 2), 
#                           labels=c("Yes","No"))
# bpq <- bpq_all[, c("SEQN", "bpq_num")]
# 
# # Merge in hypertension data
# demo_smo_hh_paq_bpq <- demo_smo_hh_paq %>% 
#   left_join(bpq, by = "SEQN") %>%
#   left_join(bp_rdg, by = "SEQN") %>%
#   left_join(bp_meds, by = "SEQN")
# 
# bp_med15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BPQ_I.XPT', 
#                      col_select = c(SEQN, BPQ050A))
# bp_med17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPQ_J.XPT', 
#                      col_select = c(SEQN, BPQ050A))
# bp_meds <- rbind(bp_med15, bp_med17)
# # Nested questionnaire so those with missing converted to No's
# bp_meds$BPQ050A[is.na(bp_meds$BPQ050A)] <- 0



# ###Cholesterol Readings 
# 
# #HDL Cholesterol is the good cholesterol - dropped. 
# 
# ##LDL cholesterol and triglyceride readings
# 
# ldltr_11  <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/TRIGLY_G.XPT', col_select = c(SEQN, LBDLDL, LBXTR))
# ldltr_13  <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/TRIGLY_H.XPT', col_select = c(SEQN, LBDLDL, LBXTR))
# ldltr_15  <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/TRIGLY_I.XPT', col_select = c(SEQN, LBDLDL, LBXTR))
# ldltr.all <- rbind(ldltr_11,ldltr_13,ldltr_15)
# 
# 
# ##Total Cholesterol readings
# tchol_11 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/TCHOL_G.XPT",col_select = c(SEQN,LBXTC))
# tchol_13 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/TCHOL_H.XPT",col_select = c(SEQN,LBXTC))
# tchol_15 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/TCHOL_I.XPT",col_select = c(SEQN,LBXTC))
# tchol.all <-rbind(tchol_11, tchol_13, tchol_15)
# chol.all <- merge(tchol.all,ldltr.all, by="SEQN")
# chol.all$chol_num <- ifelse((chol.all$LBDLDL>=100 | chol.all$LBXTC>200),1,2)
# chol.all$chol_num <- factor(chol.all$chol_num,c(1,2),labels=c("High","Normal"))
# chol <- chol.all[,c("SEQN","chol_num")] 
# 
# #merge to giant dataset
# demo.smo.hh.paq.bpq.bmi.chol <- merge(x=demo.smo.hh.paq.bpq.bmi,y=chol,all.x=TRUE,all.y=FALSE)
# 
# ###Glucose Readings 
# glu_11 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/GLU_G.XPT", col_select = c(SEQN, LBXGLU))
# glu_13 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/GLU_H.XPT", col_select = c(SEQN, LBXGLU))
# glu_15 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/GLU_I.XPT", col_select = c(SEQN, LBXGLU))
# glu.all <- rbind(glu_11, glu_13, glu_15)
# ##Assigned Categories for BMI and dropped NA values
# #Diabetic Categories: Normal (>100 mg/dl), Pre-diabetic (100-125 mg/dl), Diabetic (>125)
# glu.all$Diabetes <- ifelse(glu.all$LBXGLU<100,1,ifelse(glu.all$LBXGLU>=100 & glu.all$LBXGLU < 126,2,3))
# glu.all$Diabetes[is.na(glu.all$LBXGLU)] <- 4
# glu.all$Diabetes <- factor(glu.all$Diabetes,levels=c(1:4),labels=c("Normal","Prediabetic","Diabetic","Missing"))
# 
# demo.smo.hh.bpq.bmi.chol.glu <- merge(x=demo.smo.hh.paq.bpq.bmi.chol,y=glu.all,all.x=TRUE,all.y=FALSE)

# ## BMI Readings 
# bmx_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BMX_I.XPT',
#                    col_select = c(SEQN, BMXBMI))
# bmx_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT',
#                    col_select = c(SEQN, BMXBMI))
# bmx_all <- rbind(bmx_15, bmx_17)
# ##Assigned Categories for BMI and dropped NA values
# #Underweight (<18.5), Normal (18.5 - <25), Overweight (25 - <30), Obese (>=30) 
# bmx_all$BMI <- ifelse(bmx_all$BMXBMI < 18.5, 1, 
#                       ifelse(bmx_all$BMXBMI >= 18.5 & bmx_all$BMXBMI < 25, 2,
#                              ifelse(bmx_all$BMXBMI >= 25 & bmx_all$BMXBMI < 30, 3, 4)))
# bmx_all$BMI<- factor(bmx_all$BMI, levels = c(1:4), 
#                      labels=c("Underweight", "Normal", "Overweight", "Obese"))
# demo_smo_hh_paq_bpq_bmi <- merge(x = demo_smo_hh_paq_bpq, y = bmx_all,
#                                  all.x = TRUE, all.y = FALSE)

# #==================== HEALTHY EATING INDEX =====================================
# 
# heiadult_nhanes <- read.sas7bdat('/Users/brianajoy/OneDrive\ -\ Harvard\ University/Migrated-P-Drive/NHANES/hei_avg.sas7bdat')
# heiadultvars <- heiadult_nhanes[, c(1, 18)]
# 
# cvdrisk_all <- merge(x = demo_smo_hh_paq_bpq, y = heiadultvars, by="SEQN", 
#                      all.x = TRUE, all.y = FALSE)
# cvdrisk_all$FED_ASSIST[is.na(cvdrisk_all$FED_ASSIST)] <- 2 #turn missing into no's
# 
# #Create HEI Quantile categorical variable
# hei <- cvdrisk_all$HEI2015_TOTAL_SCORE
# hei_q1 <- quantile(hei, 0.25, na.rm=TRUE)
# hei_q2 <- quantile(hei, 0.5, na.rm=TRUE)
# hei_q3 <- quantile(hei, 0.75, na.rm=TRUE)
# cvdrisk_all$hei_cat <- ifelse(hei <= hei_q1, 1,
#                               ifelse(hei > hei_q1 & hei < hei_q2, 2,
#                                      ifelse(hei >= hei_q2 & hei < hei_q3, 3, 4)))



# #Use the Survey design to create your table 1's
# 
# nhanes.cvdpool <- svydesign(id = ~SDMVPSU,
#                             weights = ~WTMEC6YR,
#                             strata=~SDMVSTRA,
#                             nest = TRUE,
#                             data = cvdrisk.all)
# 
# 
# #Example running for each of the CVD risk factors 
# svymean(~smoke_num+Household_Smoker+bpq_num+physactive_num+BMI+chol_num+Diabetes , nhanes.cvdpool, na.rm=TRUE)
# 
# 
# ######## Use this to run your Latent Class Analysis
# lca.summer21 <- cvdrisk.all %>% dplyr::select(SEQN, smoke_num, Household_Smoker, physactive_num, bpq_num, BMI, hei_cat, chol_num, Diabetes)
# 
# # define function
# eff<-with(lca.summer21, cbind(smoke_num, Household_Smoker, physactive_num, bpq_num, BMI, hei_cat, chol_num, Diabetes)~1) 
# 
# min_bic <- 1000000
# for(i in 2:8){
#   lc <- poLCA(eff, lca.summer21, nclass=i, maxiter=3000, 
#               tol=1e-5, na.rm=TRUE,  
#               nrep=15, verbose=TRUE, calc.se=TRUE)
#   if(lc$bic < min_bic){
#     min_bic <- lc$bic
#     LCA_best_model<-lc
#   }
# } 
# lc <- poLCA(eff, lca.summer21, nclass=2, maxiter=3000,
#             tol=1e-5, na.rm=TRUE,nrep=20,verbose=TRUE, calc.se=TRUE)
