#! /bin/awk -f
#
# Variables that must be set externally
# - field: number of the field to be analysed
# - num_inst: number of instances for each scenario
# - avg_only: print the average values only
# - smooth: if set, exclude higer and lower values from stats (default: false, unassigned)
# Assume fields are separated by commas and spaces
# Assume the second field is ordered non-decreasing

BEGIN {
	FS = "[, \t]+"
	tmp_sum = 0
	tmp_sumsq = 0
	tmp_min = 1e+06
	tmp_max = 0
	last_inst = 0
	last_scname = 1
	s = 1
	num_inv_inst = 0	# count the number of instances that are not valid (e.g. unfeasible, no result available)
	number_regexp = "^[0-9]+([.][0-9]+)?([eE][+-][0-9]+)?$"	# reg.exp. to check for valid (positive) numbers	
	# skip header line, if necessary
	getline
}

{
	if ( $2 <= last_inst ) {
		# previous scenario has ended, store stats in array
		last_inst = last_inst - num_inv_inst
		if ( last_inst == 0 ) {
			# no valid records found, print "N/A"
			par_avg[last_scname] = "N/A"
			par_min[last_scname] = "N/A"
			par_max[last_scname] = "N/A"
			par_var[last_scname] = "N/A"
		}
		else {
			# ok, valid data
			if ( smooth )
				par_avg[last_scname] = (tmp_sum - tmp_min - tmp_max) / (last_inst-2)
			else
				par_avg[last_scname] = tmp_sum / last_inst
			par_min[last_scname] = tmp_min
			par_max[last_scname] = tmp_max
			if ( smooth )
				par_var[last_scname] = (tmp_sumsq - tmp_min^2 - tmp_max^2) / (last_inst-2) - par_avg[last_scname]^2 
			else
				par_var[last_scname] = tmp_sumsq / last_inst - par_avg[last_scname]^2
		}
		scenarios[s] = last_scname
		instances[s] = last_inst
		s++
		# reset helper variables
		tmp_sum = 0
		tmp_sumsq = 0
		tmp_min = 1e+06
		tmp_max = 0
		num_inv_inst = 0
	}
	
	val = $field
	# remove commas and blanks from val
	#~ gsub(/,/, "", val)
	#~ gsub(/[[:space:]]/, "", val)
	# chech it is a number
	if ( val ~ number_regexp ) {
		val = strtonum(val)
		tmp_sum += val
		tmp_sumsq += val^2
		if ( val > tmp_max ) 
			tmp_max = val
		if ( val < tmp_min ) 
			tmp_min = val
	}
	else {
		#print "val="val
		num_inv_inst++
	}
	
	last_scname = $1
	last_inst = $2
}

END {
	# record last scenario (should be inst = num_inst)
	last_inst = last_inst - num_inv_inst
	if ( last_inst == 0 ) {
		# no valid records found, print "N/A"
		par_avg[last_scname] = "N/A"
		par_min[last_scname] = "N/A"
		par_max[last_scname] = "N/A"
		par_var[last_scname] = "N/A"
	}
	else {
		# ok, valid data
		if ( smooth )
			par_avg[last_scname] = (tmp_sum - tmp_min - tmp_max) / (last_inst-2)
		else
			par_avg[last_scname] = tmp_sum / last_inst
		par_min[last_scname] = tmp_min
		par_max[last_scname] = tmp_max
		if ( smooth )
			par_var[last_scname] = (tmp_sumsq - tmp_min^2 - tmp_max^2) / (last_inst-2) - par_avg[last_scname]^2 
		else
			par_var[last_scname] = tmp_sumsq / last_inst - par_avg[last_scname]^2
	}
	scenarios[s] = last_scname
	instances[s] = last_inst
	s++
	
	# print all data
	if ( avg_only == 1 ) {
		for (i=1; i<s; i++)
			print par_avg[scenarios[i]]
	} else {
		print "scenario, v.medio, massimo, minimo, varianza"
		for (i=1; i<s; i++)
			print scenarios[i] ", " par_avg[scenarios[i]] ", " par_max[scenarios[i]] ", " par_min[scenarios[i]] ", " par_var[scenarios[i]]
	}
}
