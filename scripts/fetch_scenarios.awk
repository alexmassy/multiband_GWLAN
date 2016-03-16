#! /bin/awk -f

BEGIN{
	scenario = 0
	n_users = 0
	n_ap = 0
	area = 0
}

{
	if($1 == "Config") { scenario = $3 }
	if($1 == "n_User") { n_users = $3 }
	if($1 == "n_AP") { n_ap = $3 }
	if($1 == "field-size_X") { 
		area = $3*$3 
		print scenario " " n_users " " n_ap " " area
	}
}
