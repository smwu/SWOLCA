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


libname fped 'P:\NHANES_WSOLCA';
libname nhanes 'P:\NHANES_WSOLCA';

/*manually import demographic data from NHANES 2015-18: demo_i, demo_j*/
/* create combined demos1518 dataset */
data fped.demos1518;
	set fped.demo_i(in=c) fped.demo_j(in=d);
	wtint4yr=wtint2yr/2;  /* divide weights by 2 b/c combining 2 cycles */
	wtmec4yr=wtmec2yr/2;
if c then cycle='2015-2016';
if d then cycle='2017-2018';
keep SEQN cycle sdmvpsu sdmvstra wtint4yr wtmec4yr wtint2yr wtmec2yr riagendr ridageyr dmdmartl dmdeduc2 indfmpir ridreth3; 

run; 
proc contents data=fped.demos1518; run;

/* Restrict to those aged 20 or above according to NHANES classification as adult */
data demo1518;
	SET fped.demos1518;
	if ridageyr>=20;
run; 
 
proc contents data=fped.fped_dravg1516; run; 
proc contents data=fped.fped_dravg1718; run; 

/*shorten food group variable names for year 17-18*/
data fped_dravg1718;
	set fped.fped_dravg1718;
	array fg {37} fg1-fg37;
	array foodgroup{37} foodgroup1--foodgroup37;
	do i=1 to 37;
		fg{i}=foodgroup{i};
	end;
	drop i foodgroup1-foodgroup37;
run;

/*shorten food group variable names for year 15-16 and reorder the columns to match 17-18*/
data fped_dravg1516;
	set fped.fped_dravg1516;
	array fg {37} fg1-fg37;
	array foodgroup{37} foodgroup1 foodgroup2 foodgroup3 foodgroup4 foodgroup5 foodgroup6
	foodgroup7 foodgroup8 foodgroup9 foodgroup10 foodgroup11 foodgroup12 foodgroup13
	foodgroup14 foodgroup15 foodgroup16 foodgroup17 foodgroup18 foodgroup19 foodgroup20
	foodgroup21 foodgroup22 foodgroup23 foodgroup24 foodgroup25 foodgroup26 foodgroup27
	foodgroup28 foodgroup29 foodgroup30 foodgroup31 foodgroup32 foodgroup33 foodgroup34
	foodgroup35 foodgroup36 foodgroup37;
	do i=1 to 37;
		fg{i}=foodgroup{i};
	end;
	drop i foodgroup1-foodgroup37;
run;

/* Combine survey cycles and drop aggregate/total consumption variables */ 
data fped1518;
	set fped_dravg1516 fped_dravg1718;
	/*drop tot_fruit tot_veg tot_red tot_starchy tot_grain tot_protein tot_meat tot_dairy;*/
	drop fg1 fg5 fg7 fg10 fg15 fg18 fg19 fg30;
run;

proc contents data=fped1518 order=varnum; run;
proc sort data=fped1518; by seqn; run;
proc sort data=demo1518; by seqn; run;  

/* Combine dietary weights pulled from NHANES website */
data dietwt_ij;
	set nhanes.dr1tot_i(keep=seqn wtdrd1 wtdr2d)
		nhanes.dr1tot_j(keep=seqn wtdrd1 wtdr2d);
run; 
proc sort data=dietwt_ij; by seqn; run; 

/* Merge demographics, food pattern components, and diet weights*/
/* Keep only those with food pattern data */
data fped_nhanesadult;
	merge demo1518 fped1518(in=fped) dietwt_ij;
	by seqn;
	if ridageyr >=20; *keep only adults over 20;
	if nrecall=1 then dietwt=wtdrd1; /*assign weights according to number of completed recalls */
	else if nrecall=2 then dietwt=wtdr2d;
	dietwt4yr = dietwt/2;  /* get combined diet weight across survey cycles */
	if fped;
run; 
proc contents data=fped_nhanesadult order=varnum; run;

/* Create categorical food columns and tertiles */
data fped_adult;
	set fped_nhanesadult;
	array fg {29} fg2--fg37;
	array bb {29} bb1-bb29;
 do i=1 to 29;
 	if fg{i}=0 then bb{i}=.;  /* 0 consumption is labeled as missing */
		else bb{i}=fg{i};
 end; 
linkvar=1;
run; 

/*IDENTIFY PERSONS WHO ARE PREGNANT OR BREASTFEEDING */

proc contents data=fped.rhq_i varnum; run; 

data fped.rhq1518;
	set fped.rhq_i fped.rhq_j;

	if rhd143=1 then preg=1; /*"Are you pregnant now?"*/
		else preg=0;
	if rhq200=1 then bfeed=1; /*"Now breastfeeding a child?"*/
		else bfeed=0;
	if preg=1 or bfeed=1 then bfpreg=1;
		else bfpreg=0;
keep seqn bfeed preg bfpreg;
run; 

proc freq data=fped.rhq1518;
	tables bfeed preg bfpreg;
run;

proc sort data=fped.rhq1518; by seqn; run; 

proc sort data= fped_adult; by seqn; run;


/*Exclusion criteria 1: Restrict to female adults */
data fped_rhq_ex1;
	merge fped.rhq1518 fped_adult(in=a); 
	by seqn; 
	if a;
	if riagendr=2;
run; 

/*Exclusion criteria 2: Restrict to living at or below 130% poverty income level */
data fped_rhq_ex2;
	set fped_rhq_ex1; 
	if indfmpir = . then delete;
	if indfmpir <=1.3 ; 
run;

/*examine breastfeed/lactating mom consumption of legumes */ 
proc sort data=fped_rhq_ex2; by bfpreg; run; 
proc surveymeans data=fped_rhq_ex2 percentile=(33,66) nobs nmiss plots=none;
	strata sdmvstra;
	cluster sdmvpsu;
	weight dietwt4yr;
	var bb10 bb22;
	by bfpreg;
	*ods output Quantiles=nhanesquant_lowF;
run; 

/*Exclusion criteria 3: Remove women who are pregnant or breastfeeding */
data fped_rhq_ex3;
	set fped_rhq_ex2;
if bfpreg=1 or bfpreg=. then delete;
run;


/*CALCULATE TERTILES*/

/*Calculate tertiles based on those not lactating or breastfeeding dataset */ 
proc surveymeans data=fped_rhq_ex3 percentile=(33,66) nobs nmiss plots=none;
	strata sdmvstra;
	cluster sdmvpsu;
	weight dietwt4yr;
	var bb1-bb29;
	ods output Quantiles=nhanesquant_lowF;
run; 

/* Create first tertile cutoff using 33rd percentile */
data wtprctl33;
	set nhanesquant_lowF;
if quantile ^= 0.33 then delete;
keep varname estimate;
run; 

/* Create second tertile cutoff using 66th percentile */
data wtprctl66;
	set nhanesquant_lowF;
if quantile ^= 0.66 then delete;
keep varname estimate;
run; 

/* Convert percentiles to wide format with a column for each food group */
proc transpose data=wtprctl33 out=widepct33 prefix=pctl33_;
    id varname;
    var estimate;
run;
proc transpose data=wtprctl66 out=widepct66 prefix=pctl66_;
    id varname;
    var estimate;
run;

data widepct33;
	set widepct33;
	linkvar=1;
run; 

data widepct66;
	set widepct66;
	linkvar=1;
run; 

/* Convert bb variables to categorical based on consumption tertiles */
data tert_nhaneslowf;
	merge fped_rhq_ex3 widepct33 widepct66;
	by linkvar;
	array fg {29} fg2--fg37;
	array pctl33_bb {29} pctl33_bb1-pctl33_bb29;
	array pctl66_bb {29} pctl66_bb1-pctl66_bb29;
	array bb {29} bb1-bb29;
 	do i=1 to 29;
 	/* 
 	If no consumption, level = 1. 
 	If consumption <= 33%, level = 2
 	If consumption in (33,66]%, level = 3
 	if consumption > 66%, level = 4.
 	*/
	 	if fg{i}=0 then bb{i}=1;
			else if fg{i} >0 and fg{i} <= pctl33_bb{i} then bb{i}=2;
				else if fg{i}> pctl33_bb{i} and fg{i} <= pctl66_bb{i} then  bb{i}=3;
					else if fg{i} > pctl66_bb{i} then bb{i}=4;
	end;
	drop pctl33_bb1-pctl33_bb29 pctl66_bb1-pctl66_bb29 fg2--fg37 i linkvar  _NAME_;
run; 
proc contents data=tert_nhaneslowf varnum; run; 


proc sort data=tert_nhanesLOWF; by seqn; run; 
proc contents data=tert_nhanesLOWF varnum; run; 

/* check distribution of race/ethnicity across cycles */
proc freq data=tert_nhanesLOWF;
	tables cycle* ridreth3/nopercent;
run; 

/* check distribution of categories for each food group */
proc surveyfreq data=tert_nhanesLOWF;
	cluster sdmvpsu;
	strata sdmvstra;
	weight dietwt4yr;
	tables bb1-bb29;
run; 

/* rename variables and save data */
data fped.tert_nhaneslowF_subset1518;
	set tert_nhanesLOWF;
	rename bb1=citrus bb2=oth_fruit bb3=fruit_juice bb4=dark_green bb5=tomatoes
	bb6=oth_red bb7=potatoes bb8=oth_starchy bb9=oth_veg bb10=leg_veg 
	bb11=whole_grain bb12=ref_grain bb13=meat bb14=cured_meats bb15=organ 
	bb16=poultry bb17=seafood_high bb18=seafood_low bb19=eggs bb20=soybean 
	bb21=nuts bb22=leg_protein bb23=milk bb24=yogurt bb25=cheese 
	bb26=oils bb27=solid_fats bb28=add_sugars bb29=drinks;
run; 

/* export data to csv file */
proc export data=fped.tert_nhaneslowF_subset1518
	outfile="nhanes1518_tert_lowF.csv"
	dbms=csv
	replace;
run;


/* OLD CODE

*relabel food group variable names for year 17-18;
data fped_dravg1718;
	set fped.fped_dravg1718;
	tot_fruit=foodgroup1;
	citrus=foodgroup2;
	oth_fruit=foodgroup3;
	fruit_juice=foodgroup4;
	tot_veg=foodgroup5;
	dark_green=foodgroup6;
	tot_red=foodgroup7;
	tomatoes=foodgroup8;
	oth_red=foodgroup9;
	tot_starchy=foodgroup10;
	potatoes=foodgroup11;
	oth_starchy=foodgroup12;
	oth_veg=foodgroup13;
	leg_veg=foodgroup14;
	tot_grain=foodgroup15;
	whole_grain=foodgroup16;
	ref_grain=foodgroup17;
	tot_protein=foodgroup18;
	tot_meat=foodgroup19;
	meat=foodgroup20;
	cured_meats=foodgroup21;
	organ=foodgroup22;
	poultry=foodgroup23;
	seafood_high=foodgroup24;
	seafood_low=foodgroup25;
	eggs=foodgroup26;
	soybean=foodgroup27;
	nuts=foodgroup28;
	leg_protein=foodgroup29;
	tot_dairy=foodgroup30;
	millk=foodgroup31;
	yogurt=foodgroup32;
	cheese=foodgroup33;
	oils=foodgroup34;
	fats=foodgroup35;
	sugars=foodgroup36;
	drinks=foodgroup37;
	drop foodgroup1-foodgroup37;
run;
proc contents data=fped_dravg1718 order=varnum; run; 


*relabel food group variable names for year 15-16;
data fped_dravg1516;
	set fped.fped_dravg1516;
	tot_fruit=foodgroup1;
	citrus=foodgroup2;
	oth_fruit=foodgroup3;
	fruit_juice=foodgroup4;
	tot_veg=foodgroup5;
	dark_green=foodgroup6;
	tot_red=foodgroup7;
	tomatoes=foodgroup8;
	oth_red=foodgroup9;
	tot_starchy=foodgroup10;
	potatoes=foodgroup11;
	oth_starchy=foodgroup12;
	oth_veg=foodgroup13;
	leg_veg=foodgroup14;
	tot_grain=foodgroup15;
	whole_grain=foodgroup16;
	ref_grain=foodgroup17;
	tot_protein=foodgroup18;
	tot_meat=foodgroup19;
	meat=foodgroup20;
	cured_meats=foodgroup21;
	organ=foodgroup22;
	poultry=foodgroup23;
	seafood_high=foodgroup24;
	seafood_low=foodgroup25;
	eggs=foodgroup26;
	soybean=foodgroup27;
	nuts=foodgroup28;
	leg_protein=foodgroup29;
	tot_dairy=foodgroup30;
	millk=foodgroup31;
	yogurt=foodgroup32;
	cheese=foodgroup33;
	oils=foodgroup34;
	fats=foodgroup35;
	sugars=foodgroup36;
	drinks=foodgroup37;
	drop foodgroup1-foodgroup37;
run;
proc contents data=fped_dravg1516 order=varnum; run; 

*/

/*
	rename foodgroup1=tot_fruit foodgroup2=citrus foodgroup3=oth_fruit foodgroup4=fruit_juice
	foodgroup5=tot_veg foodgroup6=dark_green foodgroup7=tot_red foodgroup8=tomatoes
	foodgroup9=oth_red foodgroup10=tot_starchy foodgroup11=potatoes foodgroup12=oth_starchy
	foodgroup13=oth_veg foodgroup14=leg_veg foodgroup15=tot_grain foodgroup16=whole_grain 
	foodgroup17=ref_grain foodgroup18=tot_protein foodgroup19=tot_meat foodgroup20=meat
	foodgroup21=cured_meats foodgroup22=organ foodgroup23=poultry foodgroup24=seafood_high
	foodgroup25=seafood_low foodgroup26=eggs foodgroup27=soybean foodgroup28=nuts
	foodgroup29=leg_protein foodgroup30=tot_dairy foodgroup31=milk foodgroup32=yogurt 
	foodgroup33=cheese foodgroup34=oils foodgroup35=fats foodgroup36=sugars foodgroup37=drinks;
*/
