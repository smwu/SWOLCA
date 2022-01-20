#!/bin/bash

###############################################################################

# This script will make multiple jobs to be submitted based on the variables  #

# that the user "plugs in". Once  you have put in your variables and executed #

# this script make sure that you verify that the output is correct and then   #

# run submit.job.$OUTFILE and this will submit the jobs to the cluster        #

###############################################################################

###############################################################################

# If you make changes to this file you can only change it within a unix shell #

# and not wordpad, word or any other windows text editor!!! If you do make    #

# changes with *ANY* windows text editor, this script will not work!!!        #

###############################################################################


###--- The variable "OUTFILE" will create files named what you put in the variable "OUTFILE" ---###

OUTFILE="wsOFMM"

###--- OUTFILEEXT is the extension to add to the end of your file. e.g. if you want the end of your file to be named .R then put in .R, if you want your file to be named .sas then put in .sas ---###

OUTFILEEXT=".m"

###--- MYCODE is the main code that you want to create an array with. CODE reads from the file ---###

MYCODE="wsOFMM_main.m"

CODE=`cat $MYCODE`

###--- LIB is a list of libraries that you may want to load before you set up your variables. It is a list of libraries that you want to insert in your code before you define your variables. For each library that you specify, you must add a return character after the library is defined. If you do not use this then make sure you set USELIB="0". If you do use LIB then USELIB must be set to 1 ---###

LIB=""

###--- This will tell the script to use LIB (1) or not (0) ---###

USELIB="0"

###--- These are the names of the array index variables ---###

VARNAME1="bb"
VARNAME2=""
VARNAME3=""
VARNAME4=""
VARNAME5=""
VARNAME6=""
VARNAME7=""

###--- This is the number of variables you intend to use. If you only have 5 variables then put in 5. If you do not put in this field then it may break the output file ---###

NUMVAR="1"

###--- These are the values that will be used to index the array. Each variable is seperated by a space. Do not use tabs or any other delimiting character ---###

VAR1="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"
VAR2=""
VAR3=""
VAR4=""
VAR5=""
VAR6=""
VAR7=""

### This will submit this file to the cluster to create the job array. If you want to submit to a different queue or different application, you will have to modify this variable"

BSUB="sbatch --time=1-8:00:00 matlab "



## DO NOT MODIFY ANY CODE BELOW THIS LINE OTHERWISE YOU MAY BREAK THIS SCRIPT ##

###--- Initialize arrays of the index variable values ---###

array1=( $VAR1 );
array2=( $VAR2 );
array3=( $VAR3 );
array4=( $VAR4 );
array5=( $VAR5 );
array6=( $VAR6 );
array7=( $VAR7 );

###--- Set the array length so we know how many times to loop ---###

arraylength=${#array1[@]}  # Number of elements in array1

i=1;  # Counter index
while [ "$i" -le "$arraylength" ]
do
	fromarray1=${array1["$i -1"]};  # 0-based indexing
	fromarray2=${array2["$i -1"]};
	fromarray3=${array3["$i -1"]};
	fromarray4=${array4["$i -1"]};
	fromarray5=${array5["$i -1"]};
	fromarray6=${array6["$i -1"]};
	fromarray7=${array7["$i -1"]};
	
	# Load libraries if necessary
	if [ $USELIB -eq "1" ]; then
		echo "$LIB" >> $OUTFILE$i$OUTFILEEXT
	fi
	
	# Create file $OUTFILE.$i$OUTFILEEXT with index code "$VARNAME1 = $fromarray1"
	echo "$VARNAME1 = $fromarray1" >>	$OUTFILE$i$OUTFILEEXT

	if [ $NUMVAR -ge "2" ]; then
		echo "$VARNAME2 = $fromarray2" >>	$OUTFILE$i$OUTFILEEXT
	fi
	if [ $NUMVAR -ge "3" ]; then
		echo "$VARNAME3 = $fromarray3" >>	$OUTFILE$i$OUTFILEEXT
	fi
	if [ $NUMVAR -ge "4" ]; then
		echo "$VARNAME4 = $fromarray4" >>	$OUTFILE$i$OUTFILEEXT
	fi
	if [ $NUMVAR -ge "5" ]; then
		echo "$VARNAME5 = $fromarray5" >>	$OUTFILE$i$OUTFILEEXT
	fi
	if [ $NUMVAR -ge "6" ]; then
		echo "$VARNAME6 = $fromarray6" >>	$OUTFILE$i$OUTFILEEXT
	fi
	if [ $NUMVAR -ge "7" ]; then
		echo "$VARNAME7 = $fromarray7" >>	$OUTFILE$i$OUTFILEEXT
	fi
	
	# Add command to submit indexed job file to the "submit-jobs.sh" script
	### Uncomment this code to write script that submits each job individually instead of ina job array
		# echo "$BSUB $OUTFILE$i$OUTFILEEXT" >> submit-jobs.sh	

	# Append main code after index code in $OUTFILE$i$OUTFILEEXT
	echo "$CODE" >> $OUTFILE$i$OUTFILEEXT
	echo "Creating file $OUTFILE$i$OUTFILEEXT"

	let "i = $i +1";
done

# Gives you full access to "submit-jobs.sh" script while protecting against access from other users
### Uncomment this code to change execution priveleges for the "submit-jobs.sh" script
      # chmod 700 submit-jobs.sh

