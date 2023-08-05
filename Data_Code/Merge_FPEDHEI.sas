/********************************************/
/* NHANES - FPED 							*/
/* Task: Create tertile variables 			*/
/* Expand dataset to include 2017-18		*/
/********************************************/

dm log "clear";
dm output "clear";

ods preferences;
ods html close;
ods html ;

libname pride 'P:\NHANES_WSOLCA';
libname nhanes 'P:\NHANES_WSOLCA'; 

/* View, read in, and sort HEI scores */
proc contents data=pride.HEI1516_avg varnum; run; 
proc contents data=pride.HEI1718_avg varnum; run; 

data hei1518;
	set pride.hei1516_avg pride.hei1718_avg;
run; 

proc sort data=hei1518; by seqn; run; 

/* View, read in, and sort tertiles of consumption */
proc contents data=pride.tert_nhaneslowf_subset1518 varnum; run; 

data tert1518;
	set pride.tert_nhaneslowf_subset1518; 
	drop  _FREQ_;
RUN; 

proc sort data=tert1518; by seqn; run; 

/* Merge in HEI scores */
data tert_hei1518clean;
	merge hei1518(in=inhei) tert1518(in=intert);
	by SEQN;
	if intert and inhei;
run; 


/* Data check */
proc means data=tert_hei1518clean n nmiss min mean max;
	var RIDAGEYR indfmpir;
run; 

/* Save data */
data pride.hei_tert1518;
	set tert_hei1518clean;
run; 

/* export data to csv file */
proc export data=pride.hei_tert1518
	outfile="nhanes1518_hei_tert_lowF.csv"
	dbms=csv
	replace;
run;
