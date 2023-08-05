* First, download the NHANES 2015-2016 DR2TOT file and save it to your hard drive *;
* from: https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Dietary&CycleBeginYear=2015 *;
* You may need to right-click the link to the data file and select "Save target as..." *;

** TutorialUser: update this libname to reference the directory on your hard drive where you have saved the files **;
* note that the libref points to the .xpt file, NOT the directory it is saved in (unlike a typical libname) *;
libname XP xport "C:\Users\stw000\Downloads\DR2TOT_I.xpt";

* Import the xpt file and save as a temporary SAS dataset in your work directory *;
* using a data step *;
data dr2tot_i;
  set xp.dr2tot_i;
run;

* To save a permanent SAS dataset *;
** TutorialUser: update this libname to reference a directory on your hard drive where you want to save the dataset **;
libname mydata "P:\NHANES_WSOLCA";

data mydata.DR2TOT_I;
  set xp.dr2tot_i;
run;

/* Repeat for 2017-2018 data */

libname XP xport "C:\Users\stw000\Downloads\DR2TOT_J.xpt";
data dr2tot_j;
  set xp.dr2tot_j;
run;
libname mydata "P:\NHANES_WSOLCA";
data mydata.DR2TOT_J;
  set xp.dr2tot_j;
run;
