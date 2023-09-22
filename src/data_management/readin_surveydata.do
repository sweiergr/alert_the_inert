/*
	Read in survey data sets and and save as separate Stata-dta-files for each year.

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
foreach t in 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014{
	* Create file name.
	local fname_in = `"${PATH_IN_DATA}/database"'+"`t'"+".xlsx"
	local fname_out = `"${PATH_OUT_DATA}/surveydataraw"'+"`t'"
	import excel "`fname_in'", sheet("Data") firstrow clear
	* Save as dta-file in bld/out/data folder.
	save "`fname_out'", replace
}


