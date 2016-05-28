#!/bin/bash
#
# Genera il file di configurazione a partire dallo skeleton e dal dat
# Poi lancia il generatore di istanze e quindi il solver
# Itera su tutte le configurazioni presenti nel dat

ALL_CFG_FILE="gwlan-configs.dat" #esistente
INSTANCE_INI_FILE="gwlan-input.ini"
INSTANCE_SKEL_FILE="gwlan-input.ini.skel" #esistente
RES_PATH="./configurations"
NUM_INST="4"
#TIMEOUT="3600"
#REL_MIPGAP="1e-2"

CONFIGS_TO_RUN="1 2 3 4"

# $1 = skeleton file (input)
# $2 = ini file (output)


function create_ini_file ()
{	
	echo -n > $2

	cat $1 >> $2

	let j=i-1

	cat >> $2 <<EOHD

# Automatically added by the scripting system (`basename $0`):

n_AP = ${aps[$j]}
n_User = ${users[$j]}
field-size_X = ${size[$j]}
field-size_Y = ${size[$j]}
variance_present_pos = ${varpres[$j]}
variance_future_pos = ${varfut[$j]}
quantile = ${quantile[$j]}     
sigma_rayleigh = ${sigmaray[$j]}
EOHD
}

# parse configs file
cfgnr=( `cat $ALL_CFG_FILE | grep -w "Config" | awk '{print $3}'` )
users=( `cat $ALL_CFG_FILE | grep "n_User" | awk -F '=' '{print $2}'` )
aps=( `cat $ALL_CFG_FILE | grep "n_AP" | awk -F '=' '{print $2}'` )
size=( `cat $ALL_CFG_FILE | grep "field-size_X" | awk -F '=' '{print $2}'` )
varpres=( `cat $ALL_CFG_FILE | grep "variance_present_pos" | awk -F '=' '{print $2}'` )
varfut=( `cat $ALL_CFG_FILE | grep "variance_future_pos" | awk -F '=' '{print $2}'` )
quantile=( `cat $ALL_CFG_FILE | grep "quantile" | awk -F '=' '{print $2}'` )
sigmaray=( `cat $ALL_CFG_FILE | grep "sigma_rayleigh" | awk -F '=' '{print $2}'` )


for i in ${cfgnr[@]} ; do

	mkdir -p $RES_PATH/$i
	
	found="false"
	
	# only execute the specified configurations
	if [[ -n "$CONFIGS_TO_RUN" ]] ; then
		for x in $CONFIGS_TO_RUN ; do
			if [[ "$x" -eq "$i" ]] ; then
				found="true";
				break;
			fi
		done
		if [[ "$found" == "false" ]] ; then
			continue;
		fi
	fi
	
	echo "Running configuration number: $i"
	
	count="0"

	# iterate onrequired number of instances, given current configuration (i index)
	while [[  $count -lt $NUM_INST ]]; do

		count=$[$count + 1];
		
		# create input file
		create_ini_file $INSTANCE_SKEL_FILE $INSTANCE_INI_FILE

		# launch instance maker
		./inst-maker -s $count 

		# run sim (CPLEX) program
		./gwlan_ex gwlan-input.dat gwlan-input.out 0

		# save input & output files
		mkdir -p $RES_PATH/$i/$count/
		mv gwlan-input.dat $RES_PATH/$i/$count/

		mv extra.dat $RES_PATH/$i/$count/
		mv *.out $RES_PATH/$i/$count/
   
	done # iterate on instances
  	
	mv gwlan-input.ini $RES_PATH/$i/

# iterate on configs
done 

cp $ALL_CFG_FILE $RES_PATH/

#compute statistic automatically
#cd scripts/
#./fetch_statistics.sh -p ../configurations/
#cd ..
