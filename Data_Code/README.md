# Data Code
This folder include files for cleaning and processing the NHANES data. Dietary intake and hypertension data for low-income women in the United States are pooled over two survey cycles: 2015-2016 and 2017-2018. 

### Step 1: Make food groups from 24-hour diet recall data
* Data: Dietary intake from two 24-hour recalls is summarized into Food Patterns (FP) equivalents and can be downloaded from the Food Patterns Equivalents Database (FPED) website.
  * `fped_dr1tot_1516.sas7bdat` `fped_dr1tot_1718.sas7bdat` `fped_dr2tot_1516.sas7bdat` `fped_dr2tot_1718.sas7bdat`: Total nutrient intake per day across foods and beverages for each FP component
* Code:
  * `MakeFoodGroups_nhanes.sas`: Create average intake for each individual. Run separately for each survey cycle. 
* Output:
  * `fped_dravg1516.sas7bdat` `fped_dravg1718.sas7bdat`: Average intake for each FP component across the two days, for each individual and survey cycle
 
### Step 2: Create tertiles of consumption for diet data
* Data:
  * `fped_dravg1516.sas7bdat` `fped_dravg1718.sas7bdat`: Average FP intake for individuals, from Step 1
  * `demo_i.sas7bdat` `demo_j.sas7bdat`: Demographic data from NHANES website for the two survey cycles
  * `DR1TOT_I.sas7bdat` `DR1TOT_J.sas7bdat` `dr2tot_i.sas7bdat` `dr2tot_j.sas7bdat`: Dietary sampling weights obtained from NHANES website
  * `rhq_i.sas7bdat` `rhq_j.sas7bdat`: Reproductive health data from NHANES website
* Code:
  * `fped_CreatTerts_lowF.sas`: Subset to female adults age >= 20, living at or below 185% federal poverty income level, not pregnant or breastfeeding, and with reliable dietary recall data. Create tertiles of consumption (0%, 33%, 66%, 100%) for each FP component based on consumption in the dataset. 
* Output:
  * `tert_nhaneslowf_subset1518.sas7bdat`: Demographic data, sampling weights, and categorical food consumption variables with levels none, low, medium, and high
 
### Step 3: Generate HEI scores and merge with tertiles
* Data:
  * `demo_i.sas7bdat` `demo_j.sas7bdat`: Demographic data (same as Step 2)
  * `DR1TOT_I.sas7bdat` `DR1TOT_J.sas7bdat` `dr2tot_i.sas7bdat` `dr2tot_j.sas7bdat`: Diet sampling weights (same as Step 2)
  * `fped_dr1tot_1516.sas7bdat` `fped_dr1tot_1718.sas7bdat` `fped_dr2tot_1516.sas7bdat` `fped_dr2tot_1718.sas7bdat`: Total daily nutrient intake (same as Step 1)
  * `tert_nhaneslowf_subset1518.sas7bdat`: Categorized food consumption variables, from Step 2
* Code:
  * `generateHEI.sas`: Generate Healthy Eating Index (HEI) 2015 scores for all individuals. Run separately for each survey cycle. 
  * `Merge_FPEDHEI.sas`: Merge HEI scores into the tertiles of consumption dataset
* Output:
  * `hei2015_year1516.sas7bdat` `hei2015_year1718.sas7bdat`: HEI-2015 scores for both days for all individuals
  * `hei1516_avg.sas7bdat` `hei1718_avg.sas7bdat`: HEI-2015 scores averaged across days for all individuals
  * `hei_tert1518.sas7bdat`: Final processed diet data with tertiles of consumption and HEI-2015 scores.
  * `nhanes1518_hei_tert_lowF.csv`: Final processed diet data in CSV format
 
### Step 4: Add in hypertension data
* Data:
  * `nhanes1518_hei_tert_lowF.csv`: Processed diet data, from Step 3
* Code:
  * `data_app_cleaning.R`: Define binary hypertension outcome as positive if an individual has average systolic (diastolic) blood pressure of > 130 (80) over three readings, reported hypertension, or reported hypertension medication. Obtain additional sociodemographic variables of food security, smoking, physical activity. Merge all additional variables into the dietary dataset.
* Output:
  * In the `Data` folder, `nhanes1518_adult_low_f_12jul2023.csv` contains the final processed data in CSV format with categorical diet data, binary hypertension data, and additional sociodemographic variables.
