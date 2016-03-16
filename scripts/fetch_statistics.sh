#!/bin/bash
#
# Recupera i dati delle simulazioni e produce una tabella.

NUM_PARAMS=8
PARAM[1]="multib_power"
PARAM[2]="multib_time"
PARAM[3]="nr_power"
PARAM[4]="nr_time"
PARAM[5]="opt_f_p"
PARAM[6]="multib_f_p"
PARAM[7]="nr_f_p"
PARAM[8]="max_power"

PNAME[1]="Multiband Power Consumption"
PNAME[2]="Multiband Solution CPU-Time"
PNAME[3]="Nominal Power Consumption"
PNAME[4]="Nominal Solution CPU-Time"
PNAME[5]="Optimal Future Power Consumption"
PNAME[6]="Multiband Future Power Consumption"
PNAME[7]="Nominal Future Power Consumption"
PNAME[8]="Maximum Power Consumption"

MODEL[1]=""

params_to_parse=( `seq $NUM_PARAMS` )
#~ models_to_parse=( 1 4 8 9 10 11 12 13 ) 	# netscicom/itc
models_to_parse=( 1 )			# gcdwcn
#~ plefiles_to_write=( 1 2 3 4 )

SUMMARY_FILE="summary.txt"
STATS_PREFIX="stats"

# extension for statistic files
EXT="txt"

# path to the simulation folder, default "."
PTSF="."
# perform parsing and stat computation, default "yes"
DO_PARSE="yes"
DO_COMPUTE="yes"
DO_HOPS="no"

while [[ -n "$1" ]] ; do
	case "$1" in
	"-p" )
		PTSF=$2
		shift 2
		;;
	"-i" )
		NUM_INST=$2
		shift 2
		;;
	"-c" )
		DO_PARSE="no"
		shift
		;;
	"-s" )
		DO_COMPUTE="no"
		shift
		;;
	"-h" )
		DO_COMPUTE="no"
		DO_PARSE="no"
		DO_HOPS="yes"
		shift
		;;
	* )
		echo "Usage: $0 [--help] [-p <path_to_result_folder>] [-i <number_of_instances>] [-c|-s|-h]"
		echo -e "\t -c : compute stats only, no parsing"
		echo -e "\t -s : parse file only, no stats computation"
		echo -e "\t -h : analyse hop, util, and delay numbers (implies no parse and no computation)"
		exit
		;;
	esac
done

STAT_PATH="$PTSF/$SUMMARY_FILE"
STAT_PATH_TMP="$STAT_PATH"".tmp"


# check number of instances
NUM_INST=(`ls $PTSF/1 -l | grep "d" | wc -l`)
echo "Number of instances for each configuration is " $NUM_INST


SCENARI=""
LS=`ls -v $PTSF`
for s in $LS ; do
	if [[ -d "$PTSF/$s" ]] ; then
		SCENARI+="$s "
	fi
done

# reg.exp. to check for all valid numbers
#~ number_regexp='^[0-9]+([.][0-9]+)?$'		#positive only, non exp. notation
number_regexp='^-?[0-9]+([.][0-9]+)?([eE][+-][0-9]+)?$'

if [[ "$DO_PARSE" == "yes" ]] ; then
 
echo ""
#echo "Parsing files in folder '$PTSF'"
#echo ""
#echo "Summary will be written to file '$STAT_PATH'"
#echo ""
#echo "Models to be analysed: ${models_to_parse[@]}"
#echo ""

#if [[ -e "$STAT_PATH" ]] ; then
#	echo "'$STAT_PATH' already esists!"
#	echo -n "Overwrite? [y/n] "
#	read OW
#	OW=${OW:0:1}
#	if [[ "$OW" != "y" && "$OW" != "Y" ]] ; then
#		exit
#	fi
#	echo ""
#fi

rm -f $STAT_PATH
touch $STAT_PATH
rm -f $STAT_PATH_TMP
touch $STAT_PATH_TMP

echo -n "Computing statistics..."
echo ""
for s in $SCENARI ; do
	echo -n "Parsing scenario n. $s "
	ISTANZE=`ls -F $PTSF/$s | grep / | cut -d '/' -f 1`
	for i in $ISTANZE ; do
		echo -n "."
		APP="$s, $i, "
		for m in ${models_to_parse[@]} ; do
			#~ OUTFILE="$PTSF/$s/$i/${MODEL[m]}-input.out"
			OUTFILE="$PTSF/$s/$i/gwlan-input.out"
			for w in `seq $NUM_PARAMS`; do
				VAR=`grep "^${PNAME[w]}" $OUTFILE 2>/dev/null | cut -d '=' -f 2 | sed s/W//g | sed s/s//g | tr -d "\r\n"`
				if [[ -z "$VAR" ]] ; then
					VAR="N/A"
				fi
				APP+=$VAR", "
			done
			# extract solving times; check that file does contain a number!
			#~ TIMEFILE="$PTSF/$s/$i/${MODEL[m]}-input.time"
			#~ VAR=`cat $TIMEFILE 2>/dev/null`
			#~ if [[ ! $VAR =~ $number_regexp ]] ; then
				#~ VAR="N/A"
			#~ fi
			#~ APP+=$VAR", "
		done
		APP=${APP%", "}
		echo "$APP" >> $STAT_PATH_TMP
	done
	echo " done!"
done


HEADER="scenario, istanza, "
for m in ${models_to_parse[@]} ; do
for f in ${params_to_parse[@]} ; do
	HEADER+=${PARAM[f]}", "
done
done
HEADER=${HEADER%", "}

echo $HEADER > $STAT_PATH
sort -V $STAT_PATH_TMP >> $STAT_PATH
rm -f $STAT_PATH_TMP

fi	# DO_PARSE


if [[ "$DO_COMPUTE" == "yes" ]] ; then
 
echo ""


# extract instance statistics

c=0
for m in ${models_to_parse[@]} ; do
	for f in ${params_to_parse[@]} 
	do
		STATFILE=$PTSF/$STATS_PREFIX"_"${PARAM[f]}"("${MODEL[m]}").txt"
		echo "${PARAM[f]}" > $STATFILE
		let "CAMPO = 2 + c*11 + f"
		#`awk -f ./analyser.awk -v field=$CAMPO -v num_inst=$NUM_INST $STAT_PATH >> $STATFILE`
		`awk -f ./analyser.awk -v field=$CAMPO -v num_inst=$NUM_INST -v avg_only=1 $STAT_PATH >> $STATFILE`
	done
	let "c = c + 1"
done


#build feasibility counter
cat > count_feas.awk << EOHD
#! /bin/awk -f

BEGIN{
	#print "unf_rob	unf_nr"
	unf_rob = 0
	unf_nr = 0
}

{
	if(\$8 == "UNFEASIBLE," || \$8 == "N/A,"){
		unf_rob+=1
	}
	if(\$9 == "UNFEASIBLE," || \$9 == "N/A,"){
		unf_nr+=1
	}
	if(\$2 == "$NUM_INST,"){
		print \$1 " robust = " $NUM_INST-unf_rob " ; nominal = " $NUM_INST-unf_nr
		unf_nr = 0
		unf_rob = 0
	}
}
EOHD

#extract number of feasible solutions (added by Alessio Massidda) 
`awk -f ./count_feas.awk $STAT_PATH > $PTSF/feasible_solutions.txt`

#extract data about scenarios (added by Alessio Massidda)
`awk -f ./fetch_scenarios.awk $PTSF/gwlan-configs.dat > $PTSF/scenarios.txt`


# put all stat values into a single file

ALLSTATFILE=$PTSF/$STATS_PREFIX"_ALL.txt"
XYZ=$PTSF/$STATS_PREFIX"_XYZ.txt"
echo "s" > $ALLSTATFILE
for s in $SCENARI ; do
echo "$s" >> $ALLSTATFILE
done
for m in ${models_to_parse[@]} ; do
for f in ${params_to_parse[@]} ; do
	STATFILE=$PTSF/$STATS_PREFIX"_"${PARAM[f]}"("${MODEL[m]}").txt"
	paste $ALLSTATFILE $STATFILE > $XYZ
	rm $STATFILE
	mv $XYZ $ALLSTATFILE
done
done


#create final table (like the one in the paper) (added by Alessio Massidda) 
`awk '{if($1 != "s"){ print }}' $ALLSTATFILE > $PTSF/tiny_ALLstats.txt`

#from 'feasible_solutions.txt' I'm appending the rob_FS and the nr_FS and from 'scenarios.txt' I'm getting #_UT , #_AP and area (from Alessio Massidda)
`paste <(awk '{ print }' $PTSF/tiny_ALLstats.txt) <(awk '{print " " $4 " " $8}' $PTSF/feasible_solutions.txt ) <(awk '{ print " " $2 " " $3 " " $4}' $PTSF/scenarios.txt) > $PTSF/temp_stats.txt`
`rm $PTSF/tiny_ALLstats.txt $PTSF/feasible_solutions.txt $PTSF/scenarios.txt`
`awk -f ./compute_table.awk $PTSF/temp_stats.txt > $PTSF/GOOD_TABLE.txt`


echo "done!"

fi	# DO_COMPUTE


if [[ "$DO_HOPS" == "yes" ]] ; then

echo ""
echo "Analysing hop numbers..."
echo "(beware that some values in the output files might be non-meaningful)"
#~ echo "(check whether to run twice if the folder is not clean)"

for m in ${models_to_parse[@]} ; do
for f in ${plefiles_to_write[@]} ; do
	echo -n > $PTSF/${PLE_FILE[f]}.${MODEL[m]}.full.$EXT
done
done

for m in ${models_to_parse[@]} ; do
echo "Parsing model '${MODEL[m]}'..."

cntS=0
for s in $SCENARI ; do
	echo -n "Parsing scenario n. $s..."
	
	ISTANZE=`ls -F $PTSF/$s | grep / | cut -d '/' -f 1`
	cntI=0
	for i in $ISTANZE ; do
		# compute stats on single instances (avg, max, min, stdev, hist)
		# place result in i.txt files under scenario folder
		
		for f in ${plefiles_to_write[@]} ; do
			echo -n > $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT
			echo -n > $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp

			OUTFILE="$PTSF/$s/$i/${MODEL[m]}-input.out"
			if [[ $f -ne 3 ]] ; then
				grep "^${PLE_KEY[f]}" $OUTFILE 2>/dev/null | cut -d '=' -f 2 >> $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp
			else	# GW_util needs a special manipuilation
				grep "^GW_util" $OUTFILE 2>/dev/null | cut -d '=' -f 2 | sed 's/ //g' >> $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp
				sed -i "s/\r//" $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp
				sed -i "/^0$/d" $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp
			fi
			awk -f ./simple-stats.awk $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp >> $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT

			cat $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp >> $PTSF/${PLE_FILE[f]}.${MODEL[m]}.full.$EXT
			rm $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$i.$EXT.tmp
		done
		cntI=$((cntI+1))
	done

	# extract stats(estimates) from all instances of each scenario (confidence, avg, max, min, stdev)
	# place result in txt file under scenario folder

	for f in ${plefiles_to_write[@]} ; do
		echo -n > $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$EXT
		echo -n > $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
		grep "^average" $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.*.$EXT 2>/dev/null | cut -d '=' -f 2 >> $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
		awk -f ./simple-stats.awk -v pop_size=$cntI $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp >> $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$EXT
		rm $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
		rm $PTSF/$s/${PLE_FILE[f]}.${MODEL[m]}.*.$EXT
	done
	
	cntS=$((cntS+1))
	
	echo " done!"
done

# extract global stats(estimates) from all scenarios (confidence, avg, max, min, stdev)
# place result in a single txt file under sim folder

for f in ${plefiles_to_write[@]} ; do
	echo -n > $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT
	echo -n > $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
	cntF=0
	#~ grep "^average" $PTSF/*/${PLE_FILE[f]}.${MODEL[m]}.$EXT 2>/dev/null | cut -d '=' -f 2 >> $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
	# compute average+99CI for each instance and use this value for subsequent processing
	for p in `ls $PTSF/*/${PLE_FILE[f]}.${MODEL[m]}.$EXT` ; do
		a=`grep "^average" $p 2>/dev/null | cut -d '=' -f 2`
		if [[ ${PLE_DO99[f]} -eq 1 ]] ; then
			b=`grep "^99-CI" $p 2>/dev/null | cut -d '=' -f 2`
			c=`echo "scale=6; $a+$b" | bc -q 2>/dev/null`
			echo $c >> $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
		else
			echo $a >> $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
		fi
		cntF=$((cntF+1))
	done
	awk -f ./simple-stats.awk -v pop_size=$cntF $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp >> $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT
	rm $PTSF/${PLE_FILE[f]}.${MODEL[m]}.$EXT.tmp
done

echo "... done with model '${MODEL[m]}'."
done	# models_to_parse
	
# patch per tirar fuori i valori massimi assoluti
echo 
echo "Extracting global max values..."

echo -n > $PTSF/allmaxvalues.txt
for m in ${models_to_parse[@]} ; do
for f in ${plefiles_to_write[@]} ; do
	echo "Analysing" $f ${PLE_FILE[f]} $m ${MODEL[m]}
	max=`awk -f ./simple-stats.awk $PTSF/${PLE_FILE[f]}.${MODEL[m]}.full.$EXT | grep maximum | cut -d '=' -f 2`
	echo ${PLE_FILE[f]} ${MODEL[m]} $max >> $PTSF/allmaxvalues.txt
done
done

echo "...done"

fi	# DO_HOPS

echo
echo -e "Bye!\n"

exit 0

