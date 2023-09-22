/*
	Read in supplier market share in natural gas market and save as Stata-dta-file.

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
* Create file name.
local fname_in = `"${PATH_IN_DATA}/data_macro_gas_FL.xlsx"'
local fname_out = `"${PATH_OUT_DATA}/data_macro_gas_FL"'
import excel "`fname_in'", firstrow clear
*rename var
label var size_elec "Total number of HH electricity access points (yearly)"
label var size_gas "Total number of HH gas access points (yearly)"
label var mshare_2011 "2011 share in elec market"
* encode firm
encode firm, gen(firmnew)
drop firm
rename firmnew firm
label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml
*add numberlabel
numlabel , add
order month firm
* Save as dta-file in bld/out/data folder.
save "`fname_out'", replace

