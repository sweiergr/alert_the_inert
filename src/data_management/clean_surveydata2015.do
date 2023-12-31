/*
	Clean the survey data sets and and save as separate Stata-dta-files
	and append into one survey data file.
	This file clean survey 2015.

*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

* Load Stata file generated by readin-file.
use `"${PATH_IN_DATA}/surveyraw_2015.dta"', clear

rename intnr id
drop if id==.

numlabel, add

gen age=.
replace age=R7
label var age "age or respondent"
*check
sum R7
sum age

gen agecat=1 if age<=34
replace agecat=2 if age>24 & age<=34
replace agecat=3 if age>34 & age<=44
replace agecat=4 if age>44 & age<=54
replace agecat=5 if age>54 & age<=64
replace agecat=6 if age>64
label var agecat "age category of respondent"
label define agecatl 1 "18-24" 2 "25-34" 3 "35-44" 4 "45-54" 5 "55-64" 6 "65+"
label values agecat agecatl
*check  
tab agecat
tab R7_CAT

// NOTE in variable 'hometype' we regroup "detached house, townhouse and semi detached house into HOUSE 
//      to be consistent with data collection in the year 2011
// NOTE add 0 because of year 2004
gen hometype=.
label var hometype "type of house"
replace hometype=1 if P6>=2
replace hometype=2 if P6==1
label define hometypel 0 "other" 1 "house" 2 "apartment" 
label values hometype hometypel  
*check
tab hometype
tab P6

gen socialhome=.
replace socialhome=1 if P7==1
replace socialhome=2 if P7==2
replace socialhome=3 if P7==3
label var socialhome "lives in social home"
label define socialhomel 1 "yes" 2 "no" 3 "don't know"
label values socialhome socialhomel
*check
tab socialhome
tab P7

gen ownership=.
replace ownership=1 if P8==1
replace ownership=2 if P8==2
label var ownership "home ownership"
label define ownershipl 1 "owner" 2 "tenant"
label values ownership ownershipl
*check 
tab ownership
tab P8

gen income=.
replace income=1 if P9==1
replace income=2 if P9==2
replace income=3 if P9==3
replace income=4 if P9==4
replace income=5 if P9==5
replace income=6 if P9==6
replace income=7 if P9==7
replace income=8 if P9==8
replace income=9 if P9==9
replace income=10 if P9==10
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
tab P9

rename P1 size_hh
label var size_hh "size of household"
*check
tab size_hh

rename P2 kids
label var kids "number of kids in household <14 years"
*check
tab kids

gen gender=.
replace gender=1 if P3==1
replace gender=2 if P3==2
label var gender "gender of respondent"
label define genderl 1 "male" 2 "female"
label values gender genderl
*check
tab gender
tab P3

gen education=.
replace education=1 if P4==1
replace education=2 if P4==2  | P4==3 | P4==4 | P4==5 | P4==6
replace education=3 if P4==7
replace education=4 if P4==8 | P4==9
#delimit ;
label define educationl 1 "primary education" 2 "secondary education" 
	3 "higher eudcation (non-uni)" 4 "university"
;
#delimit cr
label values education educationl
*check
tab education
tab P4

gen employment=1
replace employment=2 if P11==17
replace employment=3 if P11==18
replace employment=4 if P11==19 | P11==16
replace employment=5 if P11==14 | P11==15
replace employment=6 if P11==20
label var employment "employment status(more details available)"
#delimit ;
label define employmentl 1 "working" 2 "student" 3 "housewife" 
	4 "unemployed and other inactive" 5 "retired" 6 "don't know"
;
#delimit cr
label values employment employmentl
*check
tab P11
tab employment

*HEATING DEVICES
gen heating1=.
replace heating1=1 if EI11==1 
replace heating1=2 if EI12==1 
replace heating1=3 if EI13==1
replace heating1=4 if EI14==1
replace heating1=5 if EI15==1
replace heating1=6 if EI16==1
label var heating1 "heating device 1"
#delimit ;
label define heating1l 1 "electricity (main heating)" 2 "electricity (secondary heating)"
3 "natural gas" 4 "oil" 5 "other" 6 "don't know"
;
#delimit cr
label values heating1 heating1l
*check 
tab heating1

gen supplier=.
replace supplier=98 if LW1==25 | LW1==17
replace supplier=2 if LW1==3
replace supplier=9 if LW1==8
replace supplier=51 if LW1==23
replace supplier=10 if LW1==5
replace supplier=12 if LW1==7
replace supplier=14 if LW1==9
replace supplier=17 if LW1==11
replace supplier=20 if LW1==15
replace supplier=52 if LW1==24
replace supplier=24 if LW1==16
replace supplier=25 if LW1==6
replace supplier=28 if LW1==18
replace supplier=34 if LW1==21
replace supplier=49 if LW1==22
replace supplier=99 if LW1==26
replace supplier=45 if LW1==14
label var supplier "current electricity supplier"
#delimit ;
label define supplierl 1 "antargaz" 2 "belpower" 6 "eon" 9 "ebem" 10 "ecopower" 
	12 "ecs" 14 "elegant" 17 "eneco" 20 "essent" 24 "lampiris" 25 "edf-luminus" 
	27 "nuon" 28 "octaplus" 32 "spe-luminus" 34 "wasewind"	35 "wingas"	43 "energie2030"
	45 "eni" 49 "watz" 50 "agem TSO" 51 "eandis TSO" 52 "infrax TSO" 53 "pbe TSO"
	54 "social supplier" 95 "citypower" 96 "interelecta" 97 "iverlek" 98 "other" 
	99 "don't know"
;
#delimit cr
label values supplier supplierl
*check 
tab supplier
tab LW1

gen conscious=.
replace conscious=1 if LW10==1
replace conscious=2 if LW10==2
replace conscious=3 if LW10==3
label var conscious "have you consciously chosen your supplier"
label define consciousl 1 "yes" 2 "no" 3 "don't know"
label values conscious consciousl
*check
tab LW10
tab conscious
tab supplier if conscious==.
/* everyone was asked! (in 2014 those on socialtariff were not asked!*/


gen sw_history=.
replace sw_history=1 if LW2==1
replace sw_history=2 if LW2==2
label var sw_history "have you ever changed electricity supplier"
label define sw_historyl 1 "yes" 2 "no" 3 "don't know"
label values sw_history sw_historyl
*check
tab sw_history
tab LW2

rename LW3 sw_frequency
label var sw_frequency "how often have you changed elec supplier"
label define sw_frequencyl 0 "don't know"
label values sw_frequency sw_frequencyl
*check
tab sw_frequency
tab sw_history sw_frequency

replace sw_frequency=0 if sw_frequency==. & sw_history==1
tab sw_frequency

**********************************
* create new variable for the relative sw frequency of a respondent

*1) use max of sw_frequency for those that don't know how often they have switched
gen sw_frequency_help = sw_frequency
replace sw_frequency_help=9 if sw_frequency_help==0
tab sw_frequency_help
*2) calculate mean over the entire year (excluding zeros), i.e. conditional on having switched
egen sw_freq_mean_excl=mean(sw_frequency_help)
*3) set frequency to 0 if answered "no" or "dont't know" to sw_history or if missing
*   (as they are either on a default contract or with the TSO)
replace sw_frequency_help=0 if sw_history==2 | sw_history==3 | sw_history==.
tab sw_frequency_help
// 1000 - ok, all entries have a value
*4) calculate mean over the entire year (including zeros)
egen sw_freq_mean_incl=mean(sw_frequency_help)
*5) create final variables
gen sw_relative = sw_frequency_help / sw_freq_mean_incl
label var sw_relative "relative switching frequency per year (incl. zeros)"
tab sw_relative
gen sw_relative_excl = sw_frequency_help / sw_freq_mean_excl
label var sw_relative_excl "relative switching freq per year (excl. zeros, ie conditional on having switched)"
**********************************

gen same_supplier=.
replace same_supplier=1 if LW4==1
replace same_supplier=2 if LW4==2
label var same_supplier "elec and gas supplier is same (only HH that use gas)"
label define same_supplierl 1 "yes" 2 "no" 3 "don't know"
label values same_supplier same_supplierl
*check
tab same_supplier
tab LW4
//empty cells if HH does not use gas

gen same_aware=.
replace same_aware=1 if LW5==1
replace same_aware=2 if LW5==2
replace same_aware=3 if LW5==3
label var same_aware "did you know you can choose a different supplier for elec and gas (only HH with same supplier)"
label define same_awarel 1 "yes" 2 "no" 3 "don't know"
label values same_aware same_awarel
*check
tab same_aware
tab LW5
//empty cells if HH is not with same supplier


********************************************************************************
* Code CH101 - CH110 replies (ENERGY)  - later merge with ELEC
gen no_sw_possibility_en=.
replace no_sw_possibility_en=0 if CH101==0
replace no_sw_possibility_en=1 if CH101==1
label var no_sw_possibility_en "Haven't switched ENERGY supplier bec didn't know that I can switch"
label define no_sw_possibility_enl 0 "Not true" 1 "True"
label values no_sw_possibility_en no_sw_possibility_enl
*check
tab no_sw_possibility_en

gen no_sw_how_en=.
replace no_sw_how_en=0 if CH102==0
replace no_sw_how_en=1 if CH102==1
label var no_sw_how_en "Haven't switched ENERGY supplier bec didn't know HOW to switch"
label define no_sw_how_enl 0 "Not true" 1 "True"
label values no_sw_how_en no_sw_how_enl
*check
tab no_sw_how_en

gen no_sw_cheaper_en=.
replace no_sw_cheaper_en=0 if CH103==0
replace no_sw_cheaper_en=1 if CH103==1
label var no_sw_cheaper_en "Haven't switched ENERGY supplier bec didn't receive a better offer"
label define no_sw_cheaper_enl 0 "Not true" 1 "True"
label values no_sw_cheaper_en no_sw_cheaper_enl
*check
tab no_sw_cheaper_en

gen no_sw_comparison_en=.
replace no_sw_comparison_en=0 if CH104==0
replace no_sw_comparison_en=1 if CH104==1
label var no_sw_comparison_en "Haven't switched ENERGY supplier bec difficult to compare suppliers"
label define no_sw_comparison_enl 0 "Not true" 1 "True"
label values no_sw_comparison_en no_sw_comparison_enl
*check
tab no_sw_comparison_en

gen no_sw_effort_en=.
replace no_sw_effort_en=0 if CH105==0
replace no_sw_effort_en=1 if CH105==1
label var no_sw_effort_en "Haven't switched ENERGY supplier bec switching involves too much effort"
label define no_sw_effort_enl 0 "Not true" 1 "True"
label values no_sw_effort_en no_sw_effort_enl
*check
tab no_sw_effort_en

gen no_sw_feelgood_en=.
replace no_sw_feelgood_en=0 if CH106==0
replace no_sw_feelgood_en=1 if CH106==1
label var no_sw_feelgood_en "Haven't switched ENERGY supplier bec feel well with current supplier"
label define no_sw_feelgood_enl 0 "Not true" 1 "True"
label values no_sw_feelgood_en no_sw_feelgood_enl
*check
tab no_sw_feelgood_en

gen no_sw_interest_en=.
replace no_sw_interest_en=0 if CH107==0
replace no_sw_interest_en=1 if CH107==1
label var no_sw_interest_en "Haven't switched ENERGY supplier bec not interested in a switch"
label define no_sw_interest_enl 0 "Not true" 1 "True"
label values no_sw_interest_en no_sw_interest_enl
*check
tab no_sw_interest_en

gen no_sw_wrong_en=.
replace no_sw_wrong_en=0 if CH108==0
replace no_sw_wrong_en=1 if CH108==1
label var no_sw_wrong_en "Haven't switched ENERGY supplier bec too many things could go wrong"
label define no_sw_wrong_enl 0 "Not true" 1 "True"
label values no_sw_wrong_en no_sw_wrong_enl
*check
tab no_sw_wrong_en

gen no_sw_none_en=.
replace no_sw_none_en=0 if CH109==0
replace no_sw_none_en=1 if CH109==1
label var no_sw_none_en "Haven't switched ENERGY supplier - nothing of the above mentioned applies"
label define no_sw_none_enl 0 "Not true" 1 "True"
label values no_sw_none_en no_sw_none_enl
*check
tab no_sw_none_en

gen no_sw_dontknow_en=.
replace no_sw_dontknow_en=0 if CH110==0
replace no_sw_dontknow_en=1 if CH110==1
label var no_sw_dontknow_en "Haven't switched ENERGY supplier bec: Don't know"
label define no_sw_dontknow_enl 0 "Not true" 1 "True"
label values no_sw_dontknow_en no_sw_dontknow_enl
*check
tab no_sw_dontknow_en
********************************************************************************


********************************************************************************
* Code CH201 - CH210 replies (ELEC)  - later merge with ENERGY

gen no_sw_possibility_elec=.
replace no_sw_possibility_elec=0 if CH201==0
replace no_sw_possibility_elec=1 if CH201==1
label var no_sw_possibility_elec "Haven't switched ELEC supplier bec didn't know that I can switch"
label define no_sw_possibility_elecl 0 "Not true" 1 "True"
label values no_sw_possibility_elec no_sw_possibility_elecl
*check
tab no_sw_possibility_elec

gen no_sw_how_elec=.
replace no_sw_how_elec=0 if CH202==0
replace no_sw_how_elec=1 if CH202==1
label var no_sw_how_elec "Haven't switched ELEC supplier bec didn't know HOW to switch"
label define no_sw_how_elecl 0 "Not true" 1 "True"
label values no_sw_how_elec no_sw_how_elecl
*check
tab no_sw_how_elec

gen no_sw_cheaper_elec=.
replace no_sw_cheaper_elec=0 if CH203==0
replace no_sw_cheaper_elec=1 if CH203==1
label var no_sw_cheaper_elec "Haven't switched ELEC supplier bec didn't receive a better offer"
label define no_sw_cheaper_elecl 0 "Not true" 1 "True"
label values no_sw_cheaper_elec no_sw_cheaper_elecl
*check
tab no_sw_cheaper_elec

gen no_sw_comparison_elec=.
replace no_sw_comparison_elec=0 if CH204==0
replace no_sw_comparison_elec=1 if CH204==1
label var no_sw_comparison_elec "Haven't switched ELEC supplier bec difficult to compare suppliers"
label define no_sw_comparison_elecl 0 "Not true" 1 "True"
label values no_sw_comparison_elec no_sw_comparison_elecl
*check
tab no_sw_comparison_elec

gen no_sw_effort_elec=.
replace no_sw_effort_elec=0 if CH205==0
replace no_sw_effort_elec=1 if CH205==1
label var no_sw_effort_elec "Haven't switched ELEC supplier bec switching involves too much effort"
label define no_sw_effort_elecl 0 "Not true" 1 "True"
label values no_sw_effort_elec no_sw_effort_elecl
*check
tab no_sw_effort_elec

gen no_sw_feelgood_elec=.
replace no_sw_feelgood_elec=0 if CH206==0
replace no_sw_feelgood_elec=1 if CH206==1
label var no_sw_feelgood_elec "Haven't switched ELEC supplier bec feel well with current supplier"
label define no_sw_feelgood_elecl 0 "Not true" 1 "True"
label values no_sw_feelgood_elec no_sw_feelgood_elecl
*check
tab no_sw_feelgood_elec

gen no_sw_interest_elec=.
replace no_sw_interest_elec=0 if CH207==0
replace no_sw_interest_elec=1 if CH207==1
label var no_sw_interest_elec "Haven't switched ELEC supplier bec not interested in a switch"
label define no_sw_interest_elecl 0 "Not true" 1 "True"
label values no_sw_interest_elec no_sw_interest_elecl
*check
tab no_sw_interest_elec

gen no_sw_wrong_elec=.
replace no_sw_wrong_elec=0 if CH208==0
replace no_sw_wrong_elec=1 if CH208==1
label var no_sw_wrong_elec "Haven't switched ELEC supplier bec too many things could go wrong"
label define no_sw_wrong_elecl 0 "Not true" 1 "True"
label values no_sw_wrong_elec no_sw_wrong_elecl
*check
tab no_sw_wrong_elec

gen no_sw_none_elec=.
replace no_sw_none_elec=0 if CH209==0
replace no_sw_none_elec=1 if CH209==1
label var no_sw_none_elec "Haven't switched ELEC supplier - nothing of the above mentioned applies"
label define no_sw_none_elecl 0 "Not true" 1 "True"
label values no_sw_none_elec no_sw_none_elecl
*check
tab no_sw_none_elec

gen no_sw_dontknow_elec=.
replace no_sw_dontknow_elec=0 if CH210==0
replace no_sw_dontknow_elec=1 if CH210==1
label var no_sw_dontknow_elec "Haven't switched ELEC supplier bec: Don't know"
label define no_sw_dontknow_elecl 0 "Not true" 1 "True"
label values no_sw_dontknow_elec no_sw_dontknow_elecl
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
#delimit;
replace no_sw_how=1 if no_sw_how_elec==1 | no_sw_how_en==1
;
#delimit cr
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

gen consumption=.
replace consumption=1 if EI5==1
replace consumption=2 if EI5==2
replace consumption=3 if EI5==3
replace consumption=4 if EI5==4
replace consumption=5 if EI5==5
replace consumption=6 if EI5==6
label var consumption "yearly electricity consumption (past year)"
#delimit;
label define consumptionl 1 "<900kWh" 2 "900-2350 kWh" 3 "2350-5500 kWh" 4 "5500-13750 kWh"
	5 ">13750 kWh" 6 "don't know"
;
#delimit cr
label values consumption consumptionl
*check 
tab consumption 
tab EI5
gen bill=.
gen solar=.
replace solar=1 if EI4==1
replace solar=2 if EI4==2
label var solar "HH has solar panels"
label define solarl 1 "yes" 2 "no"
label values solar solarl
*check
tab solar
tab EI4

gen energycost=.
label var energycost "How important are energy costs within total HH costs"
replace energycost=1 if EI7==4
replace energycost=2 if EI7==3
replace energycost=3 if EI7==2
replace energycost=4 if EI7==1
replace energycost=5 if EI7==5
label define energycostl 1 "not important at all" 2 "rather not important" 3 "rather important" 4 "very important" 5 "no opinion"
label values energycost energycostl
*check
tab energycost
tab EI7

gen price_change=.
replace price_change=1 if EP3==2
replace price_change=2 if EP3==1
replace price_change=3 if EP3==3
replace price_change=4 if EP3==4
label var price_change "comparing elec price HH pays today to last year"
label define price_changel 1 "increased" 2 "decreased" 3 "no change" 4 "don't know"
label values price_change price_changel
*check
tab price_change
tab EP3


gen aware_free=.
replace aware_free=1 if EI13==1
replace aware_free=2 if EI13==0
label var aware_free "do you know that ever family in FL is entitlede to free kWh (only 50%)"
label define aware_freel 1 "yes" 2 "no"
label values aware_free aware_freel
*check
tab aware_free
tab EI13
//In 2015: ALL were asked

gen info_liberalization=.
replace info_liberalization=1 if IV1==1
replace info_liberalization=2 if IV1==2
replace info_liberalization=3 if IV1==3
label var info_liberalization "do you feel sufficiently informed about liberalized energy market (only 50%)"
label define info_liberalizationl 1 "yes" 2 "no" 3 "don't know"
label values info_liberalization info_liberalizationl
*check
tab info_liberalization
tab IV1


*How did you choose your elec supplier -> several answers!
gen choice_vtest=.
replace choice_vtest=0 if LW111==0
replace choice_vtest=1 if LW111==1
label var choice_vtest "I have chosen my supplier based on V-test"
label define choice_vtestl 0 "no" 1 "yes" 
label values choice_vtest choice_vtestl
*check
tab choice_vtest 
tab LW111

gen choice_testachat=.
replace choice_testachat=0 if LW112==0
replace choice_testachat=1 if LW112==1
label var choice_testachat "I have chosen my supplier based on Test-Achats"
label define choice_testachatl 0 "no" 1 "yes" 
label values choice_testachat choice_testachatl
*check
tab choice_testachat  
tab LW112

gen choice_comparison=.
replace choice_comparison=0 if LW113==0
replace choice_comparison=1 if LW113==1
label var choice_comparison "I have chosen my supplier based on another price_comparison"
label define choice_comparisonl 0 "no" 1 "yes" 
label values choice_comparison choice_comparisonl
*check
tab choice_comparison 
tab LW113

gen choice_group=.
replace choice_group=0 if LW114==0
replace choice_group=1 if LW114==1
label var choice_group  "I have chosen my supplier based on group purchase"
label define choice_groupl 0 "no" 1 "yes" 
label values choice_group choice_groupl
*check
tab choice_group 
tab LW114

gen choice_contacted=.
replace choice_contacted=0 if LW115==0
replace choice_contacted=1 if LW115==1
label var choice_contacted "I have chosen my supplier based on contact by supplier"
label define choice_contactedl 0 "no" 1 "yes" 
label values choice_contacted choice_contactedl
*check
tab choice_contacted 
tab LW115

gen choice_people=.
replace choice_people=0 if LW116==0
replace choice_people=1 if LW116==1
label var choice_people "I have chosen my supplier based experiences from other people"
label define choice_peoplel 0 "no" 1 "yes" 
label values choice_people choice_peoplel
*check
tab choice_people 
tab LW116

gen choice_other=.
replace choice_other=0 if LW117==0
replace choice_other=1 if LW117==1
label var choice_other "I have chosen my supplier based on other"
label define choice_otherl 0 "no" 1 "yes" 
label values choice_other choice_otherl
*check
tab choice_other 
tab LW117

*construct a price comparison variable
gen choice_pricetool =.
replace choice_pricetool=0 if choice_comparison==0 & choice_testachat==0 & choice_vtest==0
replace choice_pricetool=1 if choice_comparison==1 | choice_testachat==1 | choice_vtest==1
label var choice_pricetool "Chosen supplier based on a price comparison (not original question)"
label define choice_pricetooll 0 "no" 1 "yes" 
label values choice_pricetool choice_pricetooll
tab choice_pricetool 

gen aware_vtest=.
replace aware_vtest=1 if LW12==1
replace aware_vtest=2 if LW12==2 | LW12==3 
label var aware_vtest "do you know v-test (only if knew VREG (in 2012-13 only 50%))"
label define aware_vtestl 1 "yes" 2 "no"
label value aware_vtest aware_vtestl
*check
tab aware_vtest
tab LW12
//NOTE: in 2014-2016, this was asked to EVERYONE 
//		- in 2012 and 2013 only to those that indicated to know vreg were asked
//		  (which was asked to 50% in 2012 and 2013)
//		- in 2009-2011 only to those that indicated 
//			* to know vreg (which was asked to 100% 
//			* to have visited the website

gen done_vtest=.
replace done_vtest=1 if LW13==1
replace done_vtest=2 if LW13==2 | LW13==3 
label var done_vtest "have you done v-test (only if knew VREG (in 2012-13 only 50%) and others)"
label define done_vtestl 1 "yes" 2 "no"
label value done_vtest done_vtestl
*check
tab done_vtest
tab LW13

*adapt: some have done it (see choice_vtest but have not reported it here)
replace done_vtest=1 if choice_vtest==1 //0 additional
tab done_vtest

gen switched_vtest=.
replace switched_vtest=1 if LW15==1
replace switched_vtest=2 if LW15==2 | LW15==3
label var switched_vtest "Have you based on v-test choosen another supplier"
label define switched_vtestl 1 "yes" 2 "no"
label value switched_vtest switched_vtestl
*check
tab switched_vtest
tab LW15

*what is the difference to above
tab choice_vtest switched_vtest
tab done_vtest switched_vtest 
tab done_vtest choice_vtest

gen done_group=.
replace done_group=1 if LW16==1
replace done_group=2 if LW16==2
replace done_group=3 if LW16==3
label var done_group "have you participated in a group purchase for energy?"
label define done_groupl 1 "yes" 2 "no" 3 "don't know"
label value done_group done_groupl
*check
tab done_group
tab LW16
tab choice_group done_group

*adapt: some have done it (see choice_group) but have not reported it here
replace done_group=1 if choice_group==1
tab done_group 
 
gen switched_group=.
replace switched_group=1 if LW17==1
replace switched_group=2 if LW17==2 | LW17==3
label var switched_group "Have you based on group purchase choosen another supplier"
label define switched_groupl 1 "yes" 2 "no"
label value switched_group switched_groupl
*check
tab switched_group
tab LW17
 
*adapt: some have done it (see choice_group) but have not reported it here
replace switched_group=1 if choice_group==1
tab switched_group

gen sw_intention=.
replace sw_intention=1 if CH4==4
replace sw_intention=2 if CH4==3
replace sw_intention=3 if CH4==2
replace sw_intention=4 if CH4==1
replace sw_intention=5 if CH4==6 | CH4==5
label var sw_intention "Do you consider switching in the next 6 months (11-14) / at maturity (08-10)"
label define sw_intentionl 1 "not at all" 2 "probably not" 3 "probably" 4 "certainly" 5 "don't know (yet)" 
label value sw_intention sw_intentionl
*check
tab sw_intention
tab CH4


gen sw_saving=.
replace sw_saving=1 if CH5==1
replace sw_saving=2 if CH5==2
replace sw_saving=3 if CH5==4
replace sw_saving=4 if CH5==3
label var sw_saving "at what yearly cost saving would you switch elec supplier?"
label define sw_savingl 1 "cites a sum" 2 "any saving is fine" 3 "don't know" 4 "I don't consider switching"
label value sw_saving sw_savingl
*check
tab sw_saving
tab CH5

rename CH5_VALUE sw_saving_nb
label var sw_saving_nb "at what exact yearly cost saving would you switch elec"
*check
tab sw_saving_nb


gen green_contract=.
replace green_contract=1 if ME1==1
replace green_contract=2 if ME1==2
replace green_contract=3 if ME1==3
label var green_contract "current contract is green (only 50% asked)"
label define green_contractl 1 "yes" 2 "no" 3 "don't know"
label values green_contract green_contractl
*check
tab ME1
tab green_contract
//only 50% of ALL HH where asked
tab green_contract conscious

gen green_confidence=.
replace green_confidence=1 if ME2==1
replace green_confidence=2 if ME2==2
replace green_confidence=3 if ME2==3
label var green_confidence "are you confident that green elec is effectively green (only HH that have green contract)"
label define green_confidencel 1 "yes" 2 "no" 3 "don't know"
label values green_confidence green_confidencel
*check
tab green_confidence
tab ME2

gen green_future=.
replace green_future=1 if ME4==1
replace green_future=2 if ME4==2
replace green_future=3 if ME4==3
label var green_future "do you consider choosing a green contract in the future (only HH without green contract)"
label define green_futurel 1 "yes" 2 "no" 3 "don't know"
label values green_future green_futurel
*check
tab green_future
tab ME4

gen sw_cheaper=.
replace sw_cheaper=1 if grid_LW18_r1_LW18_1==1
replace sw_cheaper=2 if grid_LW18_r1_LW18_1==2
replace sw_cheaper=3 if grid_LW18_r1_LW18_1==3 
label var sw_cheaper "supplier chosen because cheaper (only if elec contract was signed)"
label define sw_cheaperl  1 "yes" 2 "no" 3 "don't know"
label values sw_cheaper sw_cheaperl
*check
tab sw_cheaper
tab grid_LW18_r1_LW18_1

gen sw_service=.
replace sw_service=1 if grid_LW18_r2_LW18_1==1
replace sw_service=2 if grid_LW18_r2_LW18_1==2
replace sw_service=3 if grid_LW18_r2_LW18_1==3 
label var sw_service "Have chosen supplier because better service (only if elec contract was signed)"
label define sw_servicel 1 "yes" 2 "no" 3 "don't know"	
label values sw_service sw_servicel
*check
tab sw_service
tab grid_LW18_r2_LW18_1

gen sw_reliable=.
replace sw_reliable=1 if grid_LW18_r3_LW18_1==1
replace sw_reliable=2 if grid_LW18_r3_LW18_1==2
replace sw_reliable=3 if grid_LW18_r3_LW18_1==3 
label var sw_reliable "Have chosen supplier because reliable (only if elec contract was signed)"
label define sw_reliablel 1 "yes" 2 "no" 3 "don't know"	
label values sw_reliable sw_reliablel
*check
tab sw_reliable
tab grid_LW18_r3_LW18_1

gen sw_solar=.
replace sw_solar=1 if grid_LW18_r4_LW18_1==1
replace sw_solar=2 if grid_LW18_r4_LW18_1==2
replace sw_solar=3 if grid_LW18_r4_LW18_1==3 
label var sw_solar "Have chosen supplier because HH has solar panels (only if elec contract was signed)"
label define sw_solarl 1 "yes" 2 "no" 3 "don't know"
label values sw_solar sw_solarl
*check
tab sw_solar
tab grid_LW18_r4_LW18_1

gen sw_nuclear=.
replace sw_nuclear=1 if grid_LW18_r5_LW18_1==1
replace sw_nuclear=2 if grid_LW18_r5_LW18_1==2
replace sw_nuclear=3 if grid_LW18_r5_LW18_1==3 
label var sw_nuclear "Have chosen supplier because respondent does not like nuclear (only if elec contract was signed)"
label define sw_nuclearl 1 "yes" 2 "no" 3 "don't know"
label values sw_nuclear sw_nuclearl
*check
tab sw_nuclear
tab grid_LW18_r5_LW18_1

gen sw_family=.
replace sw_family=1 if grid_LW18_r6_LW18_1==1
replace sw_family=2 if grid_LW18_r6_LW18_1==2
replace sw_family=3 if grid_LW18_r6_LW18_1==3 
label var sw_family "Have chosen supplier because same as family/acquaintances (only if elec contract was signed)"
label define sw_familyl 1 "yes" 2 "no" 3 "don't know"
label values sw_family sw_familyl
*check
tab sw_family
tab grid_LW18_r6_LW18_1

gen sw_green=.
replace sw_green=1 if grid_LW18_r7_LW18_1==1
replace sw_green=2 if grid_LW18_r7_LW18_1==2
replace sw_green=3 if grid_LW18_r7_LW18_1==3
label var sw_green "Have chosen supplier because green elec (only if elec contract was signed)"
label define sw_greenl 1 "yes" 2 "no" 3 "don't know"
label values sw_green sw_greenl
*check
tab sw_green
tab grid_LW18_r7_LW18_1

gen sw_additional=.
replace sw_additional=1 if grid_LW18_r8_LW18_1==1
replace sw_additional=2 if grid_LW18_r8_LW18_1==2
replace sw_additional=3 if grid_LW18_r8_LW18_1==3
label var sw_additional "Have chosen supplier because additional service (only if elec contract was signed)"
label define sw_additionall 1 "yes" 2 "no" 3 "don't know"
label values sw_additional sw_additionall
*check
tab sw_additional
tab grid_LW18_r8_LW18_1

#delimit;
keep id hometype socialhome ownership income size_hh kids gender education 
	employment supplier conscious sw_history sw_frequency 
	same_supplier green_contract //sw_reason_gen same_reason 
	sw_cheaper 	sw_service sw_reliable sw_solar sw_nuclear sw_family sw_green 
	sw_additional no_sw_possibility no_sw_how no_sw_cheaper 
	no_sw_comparison no_sw_effort no_sw_feelgood no_sw_interest no_sw_wrong no_sw_none 
	no_sw_dontknow 	green_confidence green_future	
	green_confidence green_future 
	//green_reason 	
	heating1 // heating2 heating3 meter 
	consumption solar energycost price_change
	//price_percent
	aware_free info_liberalization //aware_vreg 
	aware_vtest done_vtest same_aware sw_intention
	//sw_intention_default 
	sw_saving //sw_saving_nb 
	//socialtariff default
	//only 2014 -2016
	choice_pricetool choice_other choice_people choice_contacted choice_group 
	choice_comparison choice_testachat choice_vtest
	switched_vtest switched_group done_group
	age agecat bill
	sw_relative sw_relative_excl
;
#delimit cr
numlabel, add
* Create a year variable to identify surveys when appending them into one file.
gen year = 2015

*********************************
* Save cleaned 2015 survey as Stata file.
save `"${PATH_OUT_DATA}/surveydata2015_clean"', replace
****************************************************************************
