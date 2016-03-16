#! /bin/awk -f


BEGIN{
	space = "        "
	tab = 11
	print "|I|" space "|J|" space  "|I|/|J|" space "delta" space "FF_multib" space "FF_nom" space "PR_pres" space "PR_fut" space "multib_time" space "nr_time"
	ratio = 0
	delta = 0
	pr_pres = 0
	pr_fut = 0
	space_12 = " "
	space_13 = " "
	space_ratio = " "
	space_delta = "    "
	space_10 = "   "
	space_11 = "   " 
	space_prP = "   "
	space_prF = "   "
	space_3 = "     "
}

{
	#computing dynamic space
	for(i=0;i<tab-length($12);i++){
		space_12 = space_12 " "
	} 
	#computing dynamic space
	for(i=0;i<tab-length($13);i++){
		space_13 = space_13 " "
	} 

	ratio = $12/$13
	
	#computing dynamic space
	for(i=0;i<tab-length(ratio);i++){
		space_ratio = space_ratio " "
	} 
	
	delta = ($13*100)/$14
	
	#computing dynamic space
	for(i=0;i<tab-length(delta);i++){
		space_delta = space_delta " "
	}
	#computing dynamic space
	for(i=0;i<tab-length($10);i++){
		space_10 = space_10 " "
	}
	#computing dynamic space
	for(i=0;i<tab-length($11);i++){
		space_11 = space_11 " "
	}
	
	if($4 == "N/A" || $4 == 0){
		pr_pres = "N/A"	
	}
	else { pr_pres = $2/$4 }

	#computing dynamic space
	for(i=0;i<tab-length(pr_pres);i++){
		space_prP = space_prP " "
	}
	
	#computing pr_fut
	if($8 == 0 || $8 == "N/A"){
		pr_fut = "N/A"	
	}	
	else { pr_fut = $7/$8 }
	
	#computing dynamic space
	for(i=0;i<tab-length(pr_fut);i++){
		space_prF = space_prF " "
	}
	#computing dynamic space
	for(i=0;i<tab-length($3);i++){
		space_3 = space_3 " "
	}
	

	#print results
	print $12 space_12 $13 space_13 ratio space_ratio delta space_delta $10 space_10 $11 space_11 pr_pres space_prP pr_fut space_prF $3 space_3 $5
	
	#resetting spaces
	space_12 = " "
	space_13 = " "
	space_ratio = " "
	space_delta = "    "
	space_10 = "   "
	space_11 = "   " 
	space_prP = "   "
	space_prF = "   "
	space_3 = "     "
}

