/*
	Clean the survey data sets and and save as separate Stata-dta-files
	and append into one survey data file.
	This file cleans survey 2013.

*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

*********************************
* Load Stata file generated by readin-file.
use `"${PATH_OUT_DATA}/surveydataraw2013.dta"', clear

rename INTNR id
drop if id==.


// NOTE in variable 'hometype' we regroup "detached house, townhouse and semi detached house into HOUSE 
//      to be consistent with data collection in the year 2011
// NOTE add 0 because of year 2004
gen hometype=.
label var hometype "type of house"
replace hometype=1 if P8a==3 | P8a==4 | P8a==5
replace hometype=2 if P8a==2
label define hometypel 0 "other" 1 "house" 2 "apartment" 
label values hometype hometypel  
*check
tab hometype


rename Isdateensocialewoning socialhome
label var socialhome "lives in social home"
label define socialhomel 1 "yes" 2 "no" 3 "don't know"
label values socialhome socialhomel
*check
tab socialhome

rename P8b* ownership
label var ownership "home ownership"
label define ownershipl 1 "owner" 2 "tenant"
label values ownership ownershipl
*check 
tab ownership

rename P11 income
label var income "family net income"
#delimit ;
label define incomel 1 "< EUR 600" 2 "EUR 600-999.99" 3 "EUR 1,000-1,499.99"
	4 "EUR 1,500-1,999.99 " 5 "EUR 2,000-2,499.99" 6 "EUR 2,500-3,749.99"
	7 "EUR 3,750-4,999.99" 8 "EUR  5,000-6,249.99" 9 "> EUR 6,250"
	10 "I prefer not to tell"
;
#delimit cr
label values income incomel
*check
tab income	

rename P2 size_hh
label var size_hh "size of household"
*check
tab size_hh

rename P3 kids
label var kids "number of kids in household <14 years"
*check
tab kids

rename P5 gender
label var gender "gender of respondent"
label define genderl 1 "male" 2 "female"
label values gender genderl
*check
tab gender


gen education=.
label var education "highest diploma"
replace education=1 if P6==1 | P6==2 | P6==3
replace education=2 if P6==4 | P6==5 
replace education=3 if P6==6
replace education=4 if P6==7
replace education=5 if P6==8
#delimit ;
label define educationl 1 "primary education" 2 "secondary education" 
	3 "higher eudcation (non-uni)" 4 "university" 5 "don't know"
	;
#delimit cr
label values education educationl
*check
tab P6
tab education

gen employment=.
replace employment=1 if P7<=6
replace employment=2 if P7==8
replace employment=3 if P7==7
replace employment=4 if P7==10 | P7==11
replace employment=5 if P7==9
replace employment=6 if P7==12
label var employment "employment status (more details available)"
#delimit ;
label define employmentl 1 "working" 2 "student" 3 "housewife" 
	4 "unemployed and other inactive" 5 "retired" 6 "don't know"
;
#delimit cr
label values employment employmentl
*check
tab employment
tab P7

rename LW6 supplier
label var supplier "current electricity supplier"
#delimit ;
label define supplierl 1 "antargaz" 2 "belpower" 6 "eon" 9 "ebem" 10 "ecopower" 
	12 "ecs" 14 "elegant" 	17 "eneco" 20 "essent" 24 "lampiris" 25 "edf-luminus" 
	27 "nuon" 28 "octaplus" 32 "spe-luminus" 34 "wasewind"	35 "wingas"	43 "energie2030"
	45 "eni" 49 "watz" 50 "agem TSO" 51 "eandis TSO" 52 "infrax TSO" 53 "pbe TSO"
	54 "social supplier" 95 "citypower" 96 "interelecta" 97 "iverlek" 98 "other" 
	99 "don't know"
;
#delimit cr
label values supplier supplierl
*check 
tab supplier
*note
note supplier: luminus = edf-luminus + spe-luminus (what is the difference?)
note supplier: eni = eni + nuon 


rename LW5 contract
label var contract "have you signed a contract? (in 2011-13 only asked if potentially on a default)"
label define contractl 1 "yes (i.e. on a competitive trariff)" 2 "no (i.e. on default or social tariff)" 3 "don't know"
label values contract contractl
*check
tab contract
tab supplier contract

gen socialtariff=.
replace socialtariff=2 
replace socialtariff=1 if supplier==50 | supplier==51 | supplier==52 | supplier==53 | supplier==54 | supplier==97 
label var socialtariff "on a social tariff with elec (not an original question)"
label define socialtariffl 1 "yes" 2 "no"
label values socialtariff socialtariffl
*check
tab socialtariff

gen default=.
replace default=2
replace default=1 if contract==2 
replace default=3 if contract==3
label var default "on a default tariff (not an original question)"
label define defaultl 1 "yes" 2 "no" 3 "don't know"
label values default defaultl
*check
tab default
tab supplier default

rename LW7a sw_history
label var sw_history "have you ever changed electricity supplier"
label define sw_historyl 1 "yes" 2 "no" 3 "don't know"
label values sw_history sw_historyl
*check
tab sw_history
tab default sw_history


rename LW7 sw_frequency
label var sw_frequency "how often have you changed elec supplier"
label define sw_frequencyl 0 "don't know" 
label values sw_frequency sw_frequencyl
*check
tab sw_frequency

**********************************
* create new variable for the relative sw frequency of a respondent

*1) use max of sw_frequency for those that don't know how often they have switched
gen sw_frequency_help = sw_frequency
replace sw_frequency_help=6 if sw_frequency_help==0
tab sw_frequency_help
*2) calculate mean over the entire year (excluding zeros), i.e. conditional on having switched
egen sw_freq_mean_excl=mean(sw_frequency_help)
*3) set frequency to 0 if answered "no" or "dont't know" to sw_history or if missing
*   (as they are either on a default contract or with the TSO)
replace sw_frequency_help=0 if sw_history==2 | sw_history==3 | sw_history==.
tab sw_frequency_help
// 1023 - ok, all entries have a value
*4) calculate mean over the entire year (including zeros)
egen sw_freq_mean_incl=mean(sw_frequency_help)
*5) create final variables
gen sw_relative = sw_frequency_help / sw_freq_mean_incl
label var sw_relative "relative switching frequency per year (incl. zeros)"
tab sw_relative
gen sw_relative_excl = sw_frequency_help / sw_freq_mean_excl
label var sw_relative_excl "relative switching freq per year (excl. zeros, ie conditional on having switched)"
**********************************

rename LW3 sw_reason_gen
label var sw_reason_gen "reasons to look for new ENERGY supplier"
#delimit ;
label define sw_reason_genl 1 "reception of bill" 2 "media" 3 "offers/advertisment of others" 
	4 "contract expires" 5 "bad service" 6 "belgian fin and econ situation"
	7 "price evolution"	38 "group purchasing" 39 "simpler procedure"
	40 "reliability" 41 "green energy" 42 "according to consumption" 
	43 "no energy" 44 "new tenant" 45 "installation of solar panels" 
	46 "visit of a representative" 47 "moving/new connection"
	48 "I won't change" 49 "else (please specify)" 50 "don't know"
;
#delimit cr
label values sw_reason_gen sw_reason_genl
*check
tab sw_reason_gen

rename LW8b same_supplier
label var same_supplier "elec and gas supplier is same (only HH that use gas)"
label define same_supplierl 1 "yes" 2 "no" 3 "don't know"
label values same_supplier same_supplierl
*check
tab same_supplier
//empty cells if HH does not use gas

rename LW8d same_aware
label var same_aware "did you know you can choose a different supplier for elec and gas (only HH with same supplier)"
label define same_awarel 1 "yes" 2 "no" 3 "don't know"
label values same_aware same_awarel
*check
tab same_aware
//empty cells if HH does not use gas OR are not with same supplier

rename LW8eH same_reason
label var same_reason "reasons for choosing same supplier for elec and gas (only if aware that different supplier is possible)"
#delimit ;
label define same_reasonl 1 "advantageous" 2 "only one bill for elec and gas"
	3 "other (please specify)" 4 "convenient" 9 "don't know"
	10 "from the past" 11 "never taken the time to compare"
	12 "habit \ regular customer \ reliable" 13 "no other options"
	14 "I can't choose" 15 "automatic"
;
#delimit cr
label values same_reason same_reasonl
*check
tab same_reason

rename LW11 sw_cheaper
label var sw_cheaper "Have chosen supplier because cheaper (only if elec contract was signed)"
label define sw_cheaperl  1 "yes" 2 "no" 3 "don't know"
label values sw_cheaper sw_cheaperl
*check
tab sw_cheaper

rename BE sw_service
label var sw_service "Have chosen supplier because better service (only if elec contract was signed)"
label define sw_servicel 1 "yes" 2 "no" 3 "don't know"	
label values sw_service sw_servicel
*check
tab sw_service

rename BF sw_reliable
label var sw_reliable "Have chosen supplier because reliable (only if elec contract was signed)"
label define sw_reliablel 1 "yes" 2 "no" 3 "don't know"	
label values sw_reliable sw_reliablel
*check
tab sw_reliable

rename BG sw_solar
label var sw_solar "Have chosen supplier because HH has solar panels (only if elec contract was signed)"
label define sw_solarl 1 "yes" 2 "no" 3 "don't know"
label values sw_solar sw_solarl
*check
tab sw_solar

rename BH sw_nuclear
label var sw_nuclear "Have chosen supplier because respondent does not like nuclear (only if elec contract was signed)"
label define sw_nuclearl 1 "yes" 2 "no" 3 "don't know"
label values sw_nuclear sw_nuclearl
*check
tab sw_nuclear

rename BI sw_family
label var sw_family "supplier chosen because same as family/acquaintances (only if elec contract was signed)"
label define sw_familyl 1 "yes" 2 "no" 3 "don't know"
label values sw_family sw_familyl
*check
tab sw_family

rename BJ sw_green
label var sw_green "Have chosen supplier because green elec (only if elec contract was signed)"
label define sw_greenl 1 "yes" 2 "no" 3 "don't know"
label values sw_green sw_greenl
*check
tab sw_green

rename BK sw_additional
label var sw_additional "Have chosen supplier because additional service (only if elec contract was signed)"
label define sw_additionall 1 "yes" 2 "no" 3 "don't know"
label values sw_additional sw_additionall
*check
tab sw_additional

foreach x of varlist LW23_W2_a CJ CK CL CM CN CO{
	local i = `i'+1
	rename `x' no_sw_energy`i'
	label var no_sw_energy`i' "why have you not switched energy supplier"
	#delimit;
	label define no_sw_energy`i'l 1 "didn't know that I can switch" 2 "didn't know HOW to switch"
	3 "no cheaper possibility" 4 "difficult to compare suppliers" 5 "asks too much effort"
	6 "feel well with supplier" 7 "not interested in switching" 8 "too many things could go wrong"
	48 "none of that" 50 "don't know"
	;
	#delimit cr
	label values no_sw_energy`i' no_sw_energy`i'l
	*check
	tab no_sw_energy`i' 
}

********************************************************************************
* Code these replies (ENERGY) as dummies - later merge with ELEC

gen no_sw_possibility_en=0
label var no_sw_possibility_en "Haven't switched ENERGY supplier bec didn't know that I can switch"
#delimit;
replace no_sw_possibility_en=1 if no_sw_energy1==1 | no_sw_energy2==1 | no_sw_energy3==1 | 
no_sw_energy4==1 | no_sw_energy5==1 | no_sw_energy6==1 | no_sw_energy7==1
;
#delimit cr
#delimit;
replace no_sw_possibility_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==. 
;
#delimit cr
label define no_sw_possibility_enl 0 "Not true" 1 "True"
label values no_sw_possibility_en no_sw_possibility_enl
*check
tab no_sw_possibility_en


gen no_sw_how_en=0
label var no_sw_how_en "Haven't switched ENERGY supplier bec didn't know HOW to switch"
#delimit;
replace no_sw_how_en=1 if no_sw_energy1==2 | no_sw_energy2==2 | no_sw_energy3==2 | 
no_sw_energy4==2 | no_sw_energy5==2 | no_sw_energy6==2 | no_sw_energy7==2
;
#delimit cr
#delimit;
replace no_sw_how_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==.
;
#delimit cr
label define no_sw_how_enl 0 "Not true" 1 "True"
label values no_sw_how_en no_sw_how_enl
*check
tab no_sw_how_en

gen no_sw_cheaper_en=0
label var no_sw_cheaper_en "Haven't switched ENERGY supplier bec didn't receive a better offer"
#delimit;
replace no_sw_cheaper_en=1 if no_sw_energy1==3 | no_sw_energy2==3 | no_sw_energy3==3 | 
no_sw_energy4==3 | no_sw_energy5==3 | no_sw_energy6==3 | no_sw_energy7==3 
;
#delimit cr
#delimit;
replace no_sw_cheaper_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==. 
;
#delimit cr
label define no_sw_cheaper_enl 0 "Not true" 1 "True"
label values no_sw_cheaper_en no_sw_cheaper_enl
*check
tab no_sw_cheaper_en

gen no_sw_comparison_en=0
label var no_sw_comparison_en "Haven't switched ENERGY supplier bec difficult to compare suppliers"
#delimit;
replace no_sw_comparison_en=1 if no_sw_energy1==4 | no_sw_energy2==4 | no_sw_energy3==4 | 
no_sw_energy4==4 | no_sw_energy5==4 | no_sw_energy6==4 | no_sw_energy7==4 
;
#delimit cr
#delimit;
replace no_sw_comparison_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==. 
;
#delimit cr
label define no_sw_comparison_enl 0 "Not true" 1 "True"
label values no_sw_comparison_en no_sw_comparison_enl
*check
tab no_sw_comparison_en

gen no_sw_effort_en=0
label var no_sw_effort_en "Haven't switched ENERGY supplier bec switching involves too much effort"
#delimit;
replace no_sw_effort_en=1 if no_sw_energy1==5 | no_sw_energy2==5 | no_sw_energy3==5 | 
no_sw_energy4==5 | no_sw_energy5==5 | no_sw_energy6==5 | no_sw_energy7==5 
;
#delimit cr
#delimit;
replace no_sw_effort_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==. 
;
#delimit cr
label define no_sw_effort_enl 0 "Not true" 1 "True"
label values no_sw_effort_en no_sw_effort_enl
*check
tab no_sw_effort_en

gen no_sw_feelgood_en=0
label var no_sw_feelgood_en "Haven't switched ENERGY supplier bec feel well with current supplier"
#delimit;
replace no_sw_feelgood_en=1 if no_sw_energy1==6 | no_sw_energy2==6 | no_sw_energy3==6 | 
no_sw_energy4==6 | no_sw_energy5==6 | no_sw_energy6==6 | no_sw_energy7==6 
;
#delimit cr
#delimit;
replace no_sw_feelgood_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==. 
;
#delimit cr
label define no_sw_feelgood_enl 0 "Not true" 1 "True"
label values no_sw_feelgood_en no_sw_feelgood_enl
*check
tab no_sw_feelgood_en

gen no_sw_interest_en=0
label var no_sw_interest_en "Haven't switched ENERGY supplier bec not interested in a switch"
#delimit;
replace no_sw_interest_en=1 if no_sw_energy1==7 | no_sw_energy2==7 | no_sw_energy3==7 | 
no_sw_energy4==7 | no_sw_energy5==7 | no_sw_energy6==7 | no_sw_energy7==7
;
#delimit cr
#delimit;
replace no_sw_interest_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==.
;
#delimit cr
label define no_sw_interest_enl 0 "Not true" 1 "True"
label values no_sw_interest_en no_sw_interest_enl
*check
tab no_sw_interest_en

gen no_sw_wrong_en=0
label var no_sw_wrong_en "Haven't switched ENERGY supplier bec too many things could go wrong"
#delimit;
replace no_sw_wrong_en=1 if no_sw_energy1==8 | no_sw_energy2==8 | no_sw_energy3==8 | 
no_sw_energy4==8 | no_sw_energy5==8 | no_sw_energy6==8 | no_sw_energy7==8
;
#delimit cr
#delimit;
replace no_sw_wrong_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==.
;
#delimit cr
label define no_sw_wrong_enl 0 "Not true" 1 "True"
label values no_sw_wrong_en no_sw_wrong_enl
*check
tab no_sw_wrong_en


gen no_sw_none_en=0
label var no_sw_none_en "Haven't switched ENERGY supplier - nothing of the above mentioned applies"
#delimit;
replace no_sw_none_en=1 if no_sw_energy1==48 | no_sw_energy2==48 | no_sw_energy3==48 | 
no_sw_energy4==48 | no_sw_energy5==48 | no_sw_energy6==48 | no_sw_energy7==48 
;
#delimit cr
#delimit;
replace no_sw_none_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==.
;
#delimit cr
label define no_sw_none_enl 0 "Not true" 1 "True"
label value no_sw_none_en no_sw_none_enl
*check
tab no_sw_none_en

gen no_sw_dontknow_en=0
label var no_sw_dontknow_en "Haven't switched ENERGY supplier bec: Don't know"
#delimit;
replace no_sw_dontknow_en=1 if no_sw_energy1==50 | no_sw_energy2==50 | no_sw_energy3==50 | 
no_sw_energy4==50 | no_sw_energy5==50 | no_sw_energy6==50 | no_sw_energy7==50
;
#delimit cr
#delimit;
replace no_sw_dontknow_en=. if no_sw_energy1==. & no_sw_energy2==. & no_sw_energy3==. & 
no_sw_energy4==. & no_sw_energy5==. & no_sw_energy6==. & no_sw_energy7==.
;
#delimit cr
label define no_sw_dontknow_enl 0 "Not true" 1 "True"
label value no_sw_dontknow_en no_sw_dontknow_enl
*check
tab no_sw_dontknow_en
********************************************************************************



foreach x of varlist LW23_W2_b CQ CR CS CT CU CV{
	local j = `j'+1
	rename `x' no_sw_elec`j'
	label var no_sw_elec`j' "why have you not switched electricity supplier"
	#delimit;
	label define no_sw_elec`j'l 1 "didn't know that I can switch" 2 "didn't know HOW to switch"
	3 "no cheaper possibility" 4 "difficult to compare suppliers" 5 "asks too much effort"
	6 "feel well with supplier" 7 "not interested in switching" 8 "too many things could go wrong"
	48 "none of that" 50 "don't know"
	;
	#delimit cr
	label values no_sw_elec`j' no_sw_elec`j'l
	*check
	tab no_sw_elec`j'
}
********************************************************************************
* code these replies as dummies (ELEC) -  later merge with ENERGY

gen no_sw_possibility_elec=0
label var no_sw_possibility_elec "Haven't switched ELEC supplier bec didn't know that I can switch"
#delimit;
replace no_sw_possibility_elec=1 if no_sw_elec1==1 | no_sw_elec2==1 | no_sw_elec3==1 | 
no_sw_elec4==1 | no_sw_elec5==1 | no_sw_elec6==1 | no_sw_elec7==1 
;
#delimit cr
#delimit;
replace no_sw_possibility_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_possibility_elecl 0 "Not true" 1 "True"
label value no_sw_possibility_elec no_sw_possibility_elecl
*test
tab no_sw_possibility_elec

gen no_sw_how_elec=0
label var no_sw_how_elec "Haven't switched ELEC supplier bec didn't know HOW to switch"
#delimit;
replace no_sw_how_elec=1 if no_sw_elec1==2 | no_sw_elec2==2 | no_sw_elec3==2 | 
no_sw_elec4==2 | no_sw_elec5==2 | no_sw_elec6==2 | no_sw_elec7==2
;
#delimit cr
#delimit;
replace no_sw_how_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_how_elecl 0 "Not true" 1 "True"
label value no_sw_how_elec no_sw_how_elecl
*check
tab no_sw_how_elec

gen no_sw_cheaper_elec=0
label var no_sw_cheaper_elec "Haven't switched ELEC supplier bec didn't receive a better offer"
#delimit;
replace no_sw_cheaper_elec=1 if no_sw_elec1==3 | no_sw_elec2==3 | no_sw_elec3==3 | 
no_sw_elec4==3 | no_sw_elec5==3 | no_sw_elec6==3 | no_sw_elec7==3
;
#delimit cr
#delimit;
replace no_sw_cheaper_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_cheaper_elecl 0 "Not true" 1 "True"
label value no_sw_cheaper_elec no_sw_cheaper_elecl
*check
tab no_sw_cheaper_elec

gen no_sw_comparison_elec=0
label var no_sw_comparison_elec "Haven't switched ELEC supplier bec difficult to compare suppliers"
#delimit;
replace no_sw_comparison_elec=1 if no_sw_elec1==4 | no_sw_elec2==4 | no_sw_elec3==4 | 
no_sw_elec4==4 | no_sw_elec5==4 | no_sw_elec6==4 | no_sw_elec7==4 
;
#delimit cr
#delimit;
replace no_sw_comparison_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_comparison_elecl 0 "Not true" 1 "True"
label value no_sw_comparison_elec no_sw_comparison_elecl
*check
tab no_sw_comparison_elec

gen no_sw_effort_elec=0
label var no_sw_effort_elec "Haven't switched ELEC supplier bec switching involves too much effort"
#delimit;
replace no_sw_effort_elec=1 if no_sw_elec1==5 | no_sw_elec2==5 | no_sw_elec3==5 | 
no_sw_elec4==5 | no_sw_elec5==5 | no_sw_elec6==5 | no_sw_elec7==5
;
#delimit cr
#delimit;
replace no_sw_effort_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_effort_elecl 0 "Not true" 1 "True"
label value no_sw_effort_elec no_sw_effort_elecl
*check
tab no_sw_effort_elec

gen no_sw_feelgood_elec=0
label var no_sw_feelgood_elec "Haven't switched ELEC supplier bec feel well with current supplier"
#delimit;
replace no_sw_feelgood_elec=1 if no_sw_elec1==6 | no_sw_elec2==6 | no_sw_elec3==6 | 
no_sw_elec4==6 | no_sw_elec5==6 | no_sw_elec6==6 | no_sw_elec7==6
;
#delimit cr
#delimit;
replace no_sw_feelgood_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_feelgood_elecl 0 "Not true" 1 "True"
label value no_sw_feelgood_elec no_sw_feelgood_elecl
*check
tab no_sw_feelgood_elec 

gen no_sw_interest_elec=0
label var no_sw_interest_elec "Haven't switched ELEC supplier bec not interested in a switch"
#delimit;
replace no_sw_interest_elec=1 if no_sw_elec1==7 | no_sw_elec2==7 | no_sw_elec3==7 | 
no_sw_elec4==7 | no_sw_elec5==7 | no_sw_elec6==7 | no_sw_elec7==7 
;
#delimit cr
#delimit;
replace no_sw_interest_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_interest_elecl 0 "Not true" 1 "True"
label value no_sw_interest_elec no_sw_interest_elecl
*check
tab no_sw_interest_elec

gen no_sw_wrong_elec=0
label var no_sw_wrong_elec "Haven't switched ELEC supplier bec too many things could go wrong"
#delimit;
replace no_sw_wrong_elec=1 if no_sw_elec1==8 | no_sw_elec2==8 | no_sw_elec3==8 | 
no_sw_elec4==8 | no_sw_elec5==8 | no_sw_elec6==8 | no_sw_elec7==8 
;
#delimit cr
#delimit;
replace no_sw_wrong_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_wrong_elecl 0 "Not true" 1 "True"
label value no_sw_wrong_elec no_sw_wrong_elecl
*check
tab no_sw_wrong_elec

gen no_sw_none_elec=0
label var no_sw_none_elec "Haven't switched ELEC supplier - nothing of the above mentioned applies"
#delimit;
replace no_sw_none_elec=1 if no_sw_elec1==48 | no_sw_elec2==48 | no_sw_elec3==48 | 
no_sw_elec4==48 | no_sw_elec5==48 | no_sw_elec6==48 | no_sw_elec7==48
;
#delimit cr
#delimit;
replace no_sw_none_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_none_elecl 0 "Not true" 1 "True"
label value no_sw_none_elec no_sw_none_elecl
*check 
tab no_sw_none_elec

gen no_sw_dontknow_elec=0
label var no_sw_dontknow_elec "Haven't switched ELEC supplier bec: Don't know"
#delimit;
replace no_sw_dontknow_elec=1 if no_sw_elec1==50 | no_sw_elec2==50 | no_sw_elec3==50 | 
no_sw_elec4==50 | no_sw_elec5==50 | no_sw_elec6==50 | no_sw_elec7==50
;
#delimit cr
#delimit;
replace no_sw_dontknow_elec=. if no_sw_elec1==. & no_sw_elec2==. & no_sw_elec3==. & 
no_sw_elec4==. & no_sw_elec5==. & no_sw_elec6==. & no_sw_elec7==.
;
#delimit cr
label define no_sw_dontknow_elecl 0 "Not true" 1 "True"
label value no_sw_dontknow_elec no_sw_dontknow_elecl
*check
tab no_sw_dontknow_elec
********************************************************************************


********************************************************************************
* merge ELEC and ENERGY dummies
gen no_sw_possibility=0
label var no_sw_possibility "Haven't switched supplier bec didn't know that I can switch"
replace no_sw_possibility=1 if no_sw_possibility_elec==1 | no_sw_possibility_en==1 
replace no_sw_possibility=. if no_sw_possibility_elec==. & no_sw_possibility_en==. 
label define no_sw_possibilityl 0 "Not true" 1 "True"
label value no_sw_possibility no_sw_possibilityl
*check 
tab no_sw_possibility

gen no_sw_how=0
label var no_sw_how "Haven't switched supplier bec didn't know HOW to switch"
replace no_sw_how=1 if no_sw_how_elec==1 | no_sw_how_en==1
replace no_sw_how=. if no_sw_how_elec==. & no_sw_how_en==. 
label define no_sw_howl 0 "Not true" 1 "True"
label value no_sw_how no_sw_howl
*check
tab no_sw_how

gen no_sw_cheaper=0
label var no_sw_cheaper "Haven't switched supplier bec didn't receive a better offer"
replace no_sw_cheaper=1 if no_sw_cheaper_elec==1 | no_sw_cheaper_en==1
replace no_sw_cheaper=. if no_sw_cheaper_elec==. & no_sw_cheaper_en==. 
label define no_sw_cheaperl 0 "Not true" 1 "True"
label value no_sw_cheaper no_sw_cheaperl
*check
tab no_sw_cheaper

gen no_sw_comparison=0
label var no_sw_comparison "Haven't switched supplier bec difficult to compare suppliers"
replace no_sw_comparison=1 if no_sw_comparison_elec==1 | no_sw_comparison_en==1 
replace no_sw_comparison=. if no_sw_comparison_elec==. & no_sw_comparison_en==. 
label define no_sw_comparisonl 0 "Not true" 1 "True"
label value no_sw_comparison no_sw_comparisonl
*check
tab no_sw_comparison

gen no_sw_effort=0
label var no_sw_effort "Haven't switched supplier bec switching involves too much effort"
replace no_sw_effort=1 if no_sw_effort_elec==1 | no_sw_effort_en==1
replace no_sw_effort=. if no_sw_effort_elec==. & no_sw_effort_en==. 
label define no_sw_effortl 0 "Not true" 1 "True"
label value no_sw_effort no_sw_effortl
*check 
tab no_sw_effort

gen no_sw_feelgood=0
label var no_sw_feelgood "Haven't switched supplier bec feel well with current supplier"
replace no_sw_feelgood=1 if no_sw_feelgood_elec==1 | no_sw_feelgood_en==1
replace no_sw_feelgood=. if no_sw_feelgood_elec==. & no_sw_feelgood_en==. 
label define no_sw_feelgoodl 0 "Not true" 1 "True"
label value no_sw_feelgood no_sw_feelgoodl
*check
tab no_sw_feelgood

gen no_sw_interest=0
label var no_sw_interest "Haven't switched supplier bec not interested in a switch"
replace no_sw_interest=1 if no_sw_interest_elec==1 | no_sw_interest_en==1
replace no_sw_interest=. if no_sw_interest_elec==. & no_sw_interest_en==. 
label define no_sw_interestl 0 "Not true" 1 "True"
label value no_sw_interest no_sw_interestl
*check
tab no_sw_interest

gen no_sw_wrong=0
label var no_sw_wrong "Haven't switched supplier bec too many things could go wrong"
replace no_sw_wrong=1 if no_sw_wrong_elec==1 | no_sw_wrong_en==1
replace no_sw_wrong=. if no_sw_wrong_elec==. & no_sw_wrong_en==. 
label define no_sw_wrongl 0 "Not true" 1 "True"
label value no_sw_wrong no_sw_wrongl
*check
tab no_sw_wrong

gen no_sw_none=0
label var no_sw_none "Haven't switched supplier bec nothing of the above mentioned"
replace no_sw_none=1 if no_sw_none_elec==1 | no_sw_none_en==1
replace no_sw_none=. if no_sw_none_elec==. & no_sw_none_en==. 
label define no_sw_nonel 0 "Not true" 1 "True"
label value no_sw_none no_sw_nonel
*check
tab no_sw_none

gen no_sw_dontknow=0
label var no_sw_dontknow "Haven't switched supplier bec: Don't know"
replace no_sw_dontknow=1 if no_sw_dontknow_elec==1 | no_sw_dontknow_en==1
replace no_sw_dontknow=. if no_sw_dontknow_elec==. & no_sw_dontknow_en==. 
label define no_sw_dontknowl 0 "Not true" 1 "True"
label value no_sw_dontknow no_sw_dontknowl
*check
tab no_sw_dontknow
********************************************************************************


foreach x of varlist EI2 K L{
	local k = `k'+1
	rename `x' heating`k'
	label var heating`k' "heating device"
	#delimit ;
	label define heating`k'l 1 "electricity (main heating)" 2 "electricity (secondary heating)"
	3 "natural gas" 4 "oil" 5 "other" 6 "don't know"
	;
	#delimit cr
	label values heating`k' heating`k'l
	*check 
	tab heating`k'
}
rename EI3 meter
label var meter "separate meter night (only if electricity is main heating device)"
label define meterl 1 "yes" 2 "no" 3 "don't know"
label values meter meterl
*check
tab meter

rename EI5H consumption
label var consumption "yearly electricity consumption (past year)"
#delimit;
label define consumptionl 1 "<900kWh" 2 "900-2350 kWh" 3 "2350-5500 kWh" 4 "5500-13750 kWh"
	5 ">13750 kWh" 6 "don't know"
;
#delimit cr
label values consumption consumptionl
*check 
tab consumption 

rename EI50 solar
label var solar "HH has solar panels"
label define solarl 1 "yes" 2 "no"
label values solar solarl
*check
tab solar

rename EI7HoebelangrijkisdeQ7A energycost
label var energycost "How important are energy costs within total HH costs"
label define energycostl 1 "not important at all" 2 "rather not important" 3 "rather important" 4 "very important" 5 "no opinion"
label values energycost energycostl
*check
tab energycost

rename LW21_W2 price_change
label var price_change "comparing elec price HH pays today to last year"
label define price_changel 1 "increased" 2 "decreased" 3 "no change" 4 "don't know"
label values price_change price_changel
*check
tab price_change

rename LW21b price_percent
label var price_percent "elec price change compared to last year in % (only if change was observed)"
label define price_percentl 1 "<5%" 2 "5-10%" 3 "10-15%" 4 "15-20%" 5 ">20%" 6 "don't know"
label values price_percent price_percentl
*check
tab price_percent	

rename EI11 aware_free
label var aware_free "do you know that ever family in FL is entitlede to free kWh (only 50%)"
label define aware_freel 1 "yes" 2 "no"
label values aware_free aware_freel
*check
tab aware_free

rename IV1a info_liberalization
label var info_liberalization "do you feel sufficiently informed about liberalized energy market (only 50%)"
label define info_liberalizationl 1 "yes" 2 "no" 3 "don't know"
label values info_liberalization info_liberalizationl
*check
tab info_liberalization

rename IV3 aware_vreg
label var aware_vreg "do you know VREG (2012-13 only 50%)"
label define aware_vregl 1 "yes, I know well" 2 "yes, I know"	3 "yes, but only by name" 4 "no"
label value aware_vreg aware_vregl
*check
tab aware_vreg

rename IV8 aware_vtest
label var aware_vtest "do you know v-test (only if knew VREG (in 2012-13 only 50%))"
label define aware_vtestl 1 "yes" 2 "no"
label value aware_vtest aware_vtestl
*check
tab aware_vtest

gen done_vtest=.
replace done_vtest=1 if IV9_W2==1
replace done_vtest=2 if IV9_W2==2 | IV9_W2==3
label var done_vtest "have you done v-test (only if knew VREG (in 2012-13 only 50%) and others)"
label define done_vtestl 1 "yes" 2 "no"
label value done_vtest done_vtestl
*check
tab done_vtest
tab IV9_W2


rename LW15S sw_intention
label var sw_intention "Do you consider switching in the next 6 months (11-14) / at maturity (08-10)"
label define sw_intentionl 1 "not at all" 2 "probably not" 3 "probably" 4 "certainly" 5 "don't know (yet)"
label value sw_intention sw_intentionl
*check
tab sw_intention

rename LW24 sw_intention_default
label var sw_intention_default  "Do you consider switching in the next 6 months (11-14) / at maturity (08-10) (only if default)"
label define sw_intention_defaultl 1 "not at all" 2 "probably not" 3 "probably" 4 "certainly" 5 "don't know"
label value sw_intention_default sw_intention_defaultl
*check
tab sw_intention_default

rename LW25V sw_saving
label var sw_saving "at what yearly cost saving would you switch elec supplier?"
label define sw_savingl 1 "cites a sum" 2 "any saving is fine" 3 "don't know" 4 "I don't consider switching"
label value sw_saving sw_savingl
*check
tab sw_saving

rename LW25R sw_saving_nb
label var sw_saving_nb "at what exact yearly cost saving would you switch elec"
*check
tab sw_saving_nb

rename ME2_W2 green_contract
label var green_contract "current contract is green (only 50% asked)"
label define green_contractl 1 "yes" 2 "no" 3 "don't know"
label values green_contract green_contractl
*check
tab green_contract
tab green_contract default

rename ME3 green_confidence
label var green_confidence "are you confident that green elec is effectively green (only HH that have green contract)"
label define green_confidencel 1 "yes" 2 "no" 3 "don't know"
label values green_confidence green_confidencel
*check
tab green_confidence

rename ME4 green_future
label var green_future "do you consider choosing a green contract in the future (only HH without green contract)"
label define green_futurel 1 "yes" 2 "no" 3 "don't know"
label values green_future green_futurel
*check
tab green_future

rename ME5 green_reason
label var green_reason "why don't you consider a green contract (only HH that don't consider green contract)"
label define green_reasonl 1 "expensive" 2 "no trust in control" 3 "limited offer" 4 "no interest" 8 "other (specify)" 9 "don't know"
label values green_reason green_reasonl
*check

* Rename weight variable.
rename Gewichtformaatxxxxxx weight
* Create year variable (to facilitate merging of yearly surveys).
gen year=2013
#delimit;
keep weight id hometype socialhome ownership income size_hh kids gender education employment
	supplier contract sw_history sw_frequency same_supplier green_contract 
	sw_reason_gen same_reason sw_cheaper sw_service sw_reliable sw_solar sw_nuclear 
	sw_family sw_green sw_additional no_sw_possibility no_sw_how no_sw_cheaper
	no_sw_comparison no_sw_effort no_sw_feelgood no_sw_interest no_sw_wrong no_sw_none
	no_sw_dontknow 	green_confidence green_future green_reason heating1 heating2
	heating3 meter consumption energycost solar price_change price_percent 
	aware_free info_liberalization aware_vreg aware_vtest done_vtest same_aware
	sw_intention sw_intention_default sw_saving sw_saving_nb socialtariff default
	sw_relative sw_relative_excl
;
#delimit cr
numlabel , add

* Create a year variable to identify surveys when appending them into one file.
gen year = 2013

*********************************
* Save cleaned 2013 survey as Stata file.
save `"${PATH_OUT_DATA}/surveydata2013_clean"', replace
****************************************************************************

