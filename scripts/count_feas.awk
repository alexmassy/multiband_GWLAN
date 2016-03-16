#! /bin/awk -f

BEGIN{
	#print "unf_rob	unf_nr"
	unf_rob = 0
	unf_nr = 0
}

{
	if($8 == "UNFEASIBLE," || $8 == "N/A,"){
		unf_rob+=1
	}
	if($9 == "UNFEASIBLE," || $9 == "N/A,"){
		unf_nr+=1
	}
	if($2 == "20,"){
		print $1 " robust = " 20-unf_rob " ; nominal = " 20-unf_nr
		unf_nr = 0
		unf_rob = 0
	}
}

