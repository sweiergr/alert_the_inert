/*
	Read in Internet penetration rates for Belgium from OECD
	And dummy of regulatory activity:
	- one big awareness campaign (CREG and SPF): 2012 Sept 17-28 = 115.000 EUR) 
	- one small (online & facebook): 2015 WHEN? = 30.303 EUR)
	- (another small for SMEs: 2016 5990 EUR)
	We are only focussing on the first and largest campaign.
*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
* Import data set from Excel files in src/original_data folder.
import excel using `"${PATH_IN_DATA}/internet.xls"', sheet(extrapolated) firstrow clear
rename A month
label var fixed_broadband "historical fixed broadband penetration rates in BEL (source OECD)"
label var mobile_broadband "historical mobile broadband penetration rates in BEL (source OECD)"
label var regulator "dummy for regulator activity (+6 months)"
sort month
* Save as Stata file.
save `"${PATH_OUT_DATA}/internet"', replace
