# Supporting Libraries #
library(haven)
library(tidyverse)
library(survey)
library(psych)
library(readxl)
library(gmodels)
library(ggpubr)

## Generate NHANES data set for total adult population ##

# Demographic variables (gender, age, race, education, weight, masked variance pseudo, income poverty ratio, federal assistance)
DEMO_15 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT', 
                   col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,
                                  INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
DEMO_15 <- DEMO_15[which(DEMO_15$RIDAGEYR > 45 ), ]  #?????

DEMO_17 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT', 
                   col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,
                                  INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
DEMO_17 <- DEMO_17[which(DEMO_17$RIDAGEYR > 45 ), ]  #?????

demo.temp <- rbind(DEMO_15, DEMO_17) 

demo.temp$RIAGENDR <- factor(demo.temp$RIAGENDR, levels=c(1,2), labels=c("Male", "Female"))
race <- c("Mexican","Other Hispanic", "NH White", "NH Black", "NH Asian", "Other/Mixed")
demo.temp$RIDRETH3 <- factor(demo.temp$RIDRETH3, levels=c(1:4,6,7), labels=race)

demo.temp$EDU <- ifelse(demo.temp$DMDEDUC2 < 3,"less than HS",
                        ifelse(demo.temp$DMDEDUC2==3,"HS/GED",
                               ifelse(demo.temp$DMDEDUC2>3,
                                      "At least some college","NA")))

demo.temp$WTINT4YR <- demo.temp$WTINT2YR/2
demo.temp$WTMEC4YR <- demo.temp$WTMEC2YR/2


#PULL PARTICIPANTS ON FOOD STAMPS OR WIC
FED_15 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/FSQ_I.XPT', 
                  col_select = c(SEQN, FSQ012,FSQ162))
FED_17 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/FSQ_J.XPT', 
                  col_select = c(SEQN, FSQ012,FSQ162))

fed.all <- rbind(FED_15, FED_17) 
fed.all$FED_ASSIST <- ifelse(fed.all$FSQ012==1 | fed.all$FSQ162==1,1,2)

demo.all <- merge(x=demo.temp, y=fed.all, by = "SEQN", all.x=TRUE, all.y=FALSE)

################################
### 6 Major CVD risk factors ###
##  1. Smoking                ##
##  2. Physical Activity      ##
##  3. Hypertension           ##
##  4. Obesity                ##  
##  5. Hyper Cholesterolemia  ##
##  6. Fasting Glucose        ##                         
################################


###Smoking Questionnaires (never, former, current, second-hand)
## Smoking data
smo_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SMQ_I.XPT', 
                   col_select = c(SEQN, SMQ020,SMQ621, SMQ050Q, SMQ040))
smo_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SMQ_J.XPT', 
                   col_select = c(SEQN, SMQ020,SMQ621, SMQ050Q, SMQ040))
smo.all <- rbind(smo_15, smo_17)
# Assigned never, former, and currently smoking 2011-2016

smo.fin <- mutate(smo.all, 
                  Smoking_Status = case_when(
                    SMQ621 == 1 ~ "Never Smoked",
                    SMQ040 == 3 ~ "Former Smoker",
                    SMQ020 == 1 | SMQ621 %in% c(2:8) ~ "Current Smoker",
                    SMQ020 == 2 ~ "Light Smoker")) 
View(smo.fin)
smo.fin$smoke_num <- ifelse(smo.fin$Smoking_Status=="Never Smoked" | 
                              smo.fin$Smoking_Status=="Light Smoker", 1,
                            ifelse(smo.fin$Smoking_Status=="Former Smoker", 2, 3))
smo.fin$smoke_num <- factor(smo.fin$smoke_num, c(1:3),
                            labels=c("Not Regular","Former","Current"))

smo.keep <- smo.fin[,c("SEQN","smoke_num")]

#Merge demographic to smoking dataset
demo.smo <- merge(x=demo.all,y=smo.keep, by="SEQN",all.x=TRUE, all.y=FALSE)

## Household Smoking Data
hh_smo_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SMQFAM_I.XPT', 
                      col_select = c(SEQN, SMD470))
hh_smo_17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SMQFAM_J.XPT', 
                      col_select = c(SEQN, SMD470))

hh_smo1518 <- rbind(hh_smo_15, hh_smo_17)
hh_smo1518$Household_Smoker <- ifelse(hh_smo1518$SMD470==0,1,2)
hh_smo1518$Household_Smoker <- factor(hh_smo1518$Household_Smoker, c(1,2), 
                                      labels=c("No","Yes"))
hh_smo <- hh_smo1518[,c("SEQN","Household_Smoker")]
# Add household smoking to demo dataset
demo.smo.hh <- merge(x=demo.smo,y=hh_smo, by="SEQN",all.x=TRUE, all.y=FALSE)



###Physical activity questionnaires 
paq_11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/PAQ_G.XPT', col_select = c(SEQN, PAQ610, PAD615, PAQ640, PAD645, PAQ655, PAD660))
paq_13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/PAQ_H.XPT', col_select = c(SEQN, PAQ610, PAD615, PAQ640, PAD645, PAQ655, PAD660))
paq_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/PAQ_I.XPT', col_select = c(SEQN, PAQ610, PAD615, PAQ640, PAD645, PAQ655, PAD660))
paq.all <- rbind(paq_11, paq_13, paq_15)
#Assigned 0 Days as Not Active, 1 day as Barely Active, 2-7 Days as Active
sum1 <- function(x){
  if(all(is.na(x))){return(NA)}
  else{
    return(sum(x, na.rm = TRUE))
  }}
#NA + NA = NA , NA + n = n

paq.all <- paq.all%>%
  mutate(Min_Active1 = case_when (!is.na(PAQ610) & !is.na(PAD615) ~ PAQ610 * PAD615))
paq.all <- paq.all%>%
  mutate(Min_Active2 = case_when (!is.na(PAQ640) & !is.na(PAD645) ~ PAQ640 * PAD645))
paq.all <- paq.all%>%
  mutate(Min_Active3 = case_when (!is.na(PAQ655) & !is.na(PAD660) ~ PAQ655 * PAD660))

paq.all$Min_Active <- apply(X = paq.all[,8:10], MARGIN = 1, FUN = sum1)

paq.all <- paq.all%>%
  mutate(Phys_Active = case_when ( 
    Min_Active >= 150 ~ "Active",
    Min_Active < 150 ~ "Not Active",
  )
  )
paq.all$physactive_num <- ifelse(is.na(paq.all$Phys_Active),1,ifelse(paq.all$Phys_Active=="Not Active",2,3))
paq.all$physactive_num<- factor(paq.all$physactive_num,levels=c(1:3),labels=c("Missing","Not Active","Active"))
paq.all <- paq.all[,c("SEQN","physactive_num")]
demo.smo.hh.paq <- merge(x=demo.smo.hh,y=paq.all, by="SEQN",all.x=TRUE, all.y=FALSE)


###Blood pressure questionnaire for who had been told they had hypertension 2+ times
bpq_11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BPQ_G.XPT', col_select = c(SEQN, BPQ030))
bpq_13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BPQ_H.XPT', col_select = c(SEQN, BPQ030))
bpq_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BPQ_I.XPT', col_select = c(SEQN, BPQ030))
bpq.all <- rbind(bpq_11, bpq_13, bpq_15)

bpq.all$bpq_num <- ifelse(bpq.all$BPQ030==1,1,2)
bpq.all$bpq_num[is.na(bpq.all$bpq_num)] <- 2
bpq.all$bpq_num <- factor(bpq.all$bpq_num,levels=c(1,2), labels=c("Yes","No"))
bpq <- bpq.all[,c("SEQN","bpq_num")]
demo.smo.hh.paq.bpq <- merge(x=demo.smo.hh.paq,y=bpq,all.x=TRUE,all.y=FALSE)

### BMI Readings
bmx_11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BMX_G.XPT',col_select = c(SEQN,BMXBMI))
bmx_13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BMX_H.XPT',col_select = c(SEQN,BMXBMI))
bmx_15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BMX_I.XPT',col_select = c(SEQN,BMXBMI))
bmx.all <- rbind(bmx_11,bmx_13, bmx_15)
##Assigned Categories for BMI and dropped NA values
#Underweight (<18.5), Normal (18.5 - <25), Overweight (25 - <30), Obese (>=30) 
bmx.all$BMI <- ifelse(bmx.all$BMXBMI<18.5,1,ifelse(bmx.all$BMXBMI>=18.5 & bmx.all$BMXBMI < 25,2,ifelse(bmx.all$BMXBMI>=25 & bmx.all$BMXBMI <30,3,4)))
bmx.all$BMI<- factor(bmx.all$BMI,levels=c(1:4),labels=c("Underweight","Normal","Overweight","Obese"))
demo.smo.hh.paq.bpq.bmi <- merge(x=demo.smo.hh.paq.bpq,y=bmx.all,all.x=TRUE,all.y=FALSE)


###Cholesterol Readings 

#HDL Cholesterol is the good cholesterol - dropped. 

##LDL cholesterol and triglyceride readings

ldltr_11  <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/TRIGLY_G.XPT', col_select = c(SEQN, LBDLDL, LBXTR))
ldltr_13  <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/TRIGLY_H.XPT', col_select = c(SEQN, LBDLDL, LBXTR))
ldltr_15  <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/TRIGLY_I.XPT', col_select = c(SEQN, LBDLDL, LBXTR))
ldltr.all <- rbind(ldltr_11,ldltr_13,ldltr_15)


##Total Cholesterol readings
tchol_11 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/TCHOL_G.XPT",col_select = c(SEQN,LBXTC))
tchol_13 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/TCHOL_H.XPT",col_select = c(SEQN,LBXTC))
tchol_15 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/TCHOL_I.XPT",col_select = c(SEQN,LBXTC))
tchol.all <-rbind(tchol_11, tchol_13, tchol_15)
chol.all <- merge(tchol.all,ldltr.all, by="SEQN")
chol.all$chol_num <- ifelse((chol.all$LBDLDL>=100 | chol.all$LBXTC>200),1,2)
chol.all$chol_num <- factor(chol.all$chol_num,c(1,2),labels=c("High","Normal"))
chol <- chol.all[,c("SEQN","chol_num")] 

#merge to giant dataset
demo.smo.hh.paq.bpq.bmi.chol <- merge(x=demo.smo.hh.paq.bpq.bmi,y=chol,all.x=TRUE,all.y=FALSE)

###Glucose Readings 
glu_11 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/GLU_G.XPT", col_select = c(SEQN, LBXGLU))
glu_13 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/GLU_H.XPT", col_select = c(SEQN, LBXGLU))
glu_15 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/GLU_I.XPT", col_select = c(SEQN, LBXGLU))
glu.all <- rbind(glu_11, glu_13, glu_15)
##Assigned Categories for BMI and dropped NA values
#Diabetic Categories: Normal (>100 mg/dl), Pre-diabetic (100-125 mg/dl), Diabetic (>125)
glu.all$Diabetes <- ifelse(glu.all$LBXGLU<100,1,ifelse(glu.all$LBXGLU>=100 & glu.all$LBXGLU < 126,2,3))
glu.all$Diabetes[is.na(glu.all$LBXGLU)] <- 4
glu.all$Diabetes <- factor(glu.all$Diabetes,levels=c(1:4),labels=c("Normal","Prediabetic","Diabetic","Missing"))

demo.smo.hh.bpq.bmi.chol.glu <- merge(x=demo.smo.hh.paq.bpq.bmi.chol,y=glu.all,all.x=TRUE,all.y=FALSE)

## HEALTHY EATING INDEX 

heiadult.nhanes <- read.sas7bdat('/Users/brianajoy/OneDrive\ -\ Harvard\ University/Migrated-P-Drive/NHANES/hei_avg.sas7bdat')
heiadultvars <- heiadult.nhanes[,c(1,18)]

cvdrisk.all <- merge(x=demo.smo.hh.bpq.bmi.chol.glu,y=heiadultvars,by="SEQN",all.x=TRUE,all.y=FALSE)

#Create HEI Quantile categorical variable
hei <- cvdrisk.all$HEI2015_TOTAL_SCORE
hei.q1 <- quantile(hei,0.25,na.rm=TRUE)
hei.q2 <- quantile(hei,0.5,na.rm=TRUE)
hei.q3 <- quantile(hei,0.75,na.rm=TRUE)
cvdrisk.all$hei_cat <- ifelse(hei <= hei.q1,1,ifelse(hei>hei.q1 & hei<hei.q2,2,ifelse(hei>=hei.q2 & hei < hei.q3,3,4)))
cvdrisk.all$FED_ASSIST[is.na(cvdrisk.all$FED_ASSIST)] <- 2 #turn missing into no's

write.csv(cvdrisk.all,"nhanesadult1116_dataset14jul21.csv", row.names = FALSE)


#Use the Survey design to create your table 1's

nhanes.cvdpool <- svydesign(id = ~SDMVPSU,
                            weights = ~WTMEC6YR,
                            strata=~SDMVSTRA,
                            nest = TRUE,
                            data = cvdrisk.all)


#Example running for each of the CVD risk factors 
svymean(~smoke_num+Household_Smoker+bpq_num+physactive_num+BMI+chol_num+Diabetes , nhanes.cvdpool, na.rm=TRUE)


######## Use this to run your Latent Class Analysis
lca.summer21 <- cvdrisk.all %>% dplyr::select(SEQN, smoke_num, Household_Smoker, physactive_num, bpq_num, BMI, hei_cat, chol_num, Diabetes)

# define function
eff<-with(lca.summer21, cbind(smoke_num, Household_Smoker, physactive_num, bpq_num, BMI, hei_cat, chol_num, Diabetes)~1) #

min_bic <- 1000000
for(i in 2:8){
  lc <- poLCA(eff, lca.summer21, nclass=i, maxiter=3000, 
              tol=1e-5, na.rm=TRUE,  
              nrep=15, verbose=TRUE, calc.se=TRUE)
  if(lc$bic < min_bic){
    min_bic <- lc$bic
    LCA_best_model<-lc
  }
} 
lc <- poLCA(eff, lca.summer21, nclass=2, maxiter=3000,
            tol=1e-5, na.rm=TRUE,nrep=20,verbose=TRUE, calc.se=TRUE)


