*******************************************************************************
* This file prepares tables and graphs for the data section of the paper.
*******************************************************************************

* BASICS.
clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_ANALYSIS}/log/`1'.log"', replace
*******************************************************************************
* Aggregate data
* Load clean data from survey.
use `"${PATH_OUT_DATA}/new_master_data"', clear

*************************************
*by year display mshare and price of each of the suppliers
gen year=int(month/100)
order year
* rename the firms
decode firm, gen(newfirm)
drop firm
gen firm=substr(newfirm,4,.) 
drop newfirm
order year firm
replace firm="Eni" if firm=="ENINuon"

*************************************
*(1) by month display aggregate swrate
preserve 

*drop trash
keep month switching

*drop duplicates
quietly by month, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*prep data
sort month
gen newdate=date(string(month),"YM")
gen date=newdate
format %td date	

*graph
twoway line switching date, ///
	text(0 19237 "Campaign", color(black)) ///
	xline(18993 19359 19724 20089 20454, lcolor(red) lpattern(dash) lwidth(thin)) ///
	xline(19237, lcolor(black) lwidth(thick)) ///
	xtitle("") ///
	graphregion(color(white)) bgcolor(white) ///
	ytitle("Supplier churn rate")
	
graph export `"${PATH_OUT_FIGURES}/switching.pdf"', replace
restore
*************************************

*************************************
*(2) switching rate together with total advertisement
preserve

*drop trash
keep month firm switching ad_spending_mly

*create total advertisement
by month, sort: egen advertisement=total(ad_spending_mly)
drop firm ad_spending_mly

*drop duplicates
quietly by month advertisement, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*prep data
sort month
gen newdate=date(string(month),"YM")
gen date=newdate
format %td date	

*graph
twoway (line switching date, lwidth(thick)) (line advertisement date, lwidth(thin) yaxis(2)), ///
	text(0 19237 "Campaign", color(black)) ///
	xline(19237, lcolor(black) lwidth(thick)) ///
	xtitle("") ///
	ytitle("") ///
	graphregion(color(white)) bgcolor(white) ///
	legend(label(1 Supplier churn rate) label(2 Advertisement (in 1,000 EUR)))	

graph export `"${PATH_OUT_FIGURES}/switching_adv.pdf"', replace

corr switching advertisement

restore
*************************************

*************************************
*(3) scatter MS and advertisement
preserve

*drop trash
keep year firm contract_type mshare ad_spending_mly ad_spending_norm
drop if contract_type==1
drop contract_type

*create yearly figures
by year firm, sort: egen ad_spending_yrl=total(ad_spending_mly)
by year firm, sort: egen ad_spending_norm_yrl=total(ad_spending_norm)
by year firm, sort: egen mshare_yrl=mean(mshare)

*drop duplicates
quietly by year mshare_yrl, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

drop if firm=="Other"
*plot
separate mshare_yrl, by(firm) veryshortlabel
scatter mshare_yrl? ad_spending_yrl, ms(O ..) ///
 ytitle("Average yearly market share") ///
 xtitle("Total yearly advertisement spending") ///
 mcolor(gs0 gs4 gs7 gs10 gs14 gs16) mlcolor(black ..) msize(*1.5 ..) ///
 graphregion(color(white)) bgcolor(white) ///
 legend(pos(11) ring(0) col(1))

scatter mshare_yrl? ad_spending_yrl, xsc(log)  ms(O ..) ///
 ytitle("Average yearly market share") ///
 xtitle("Total yearly advertisement spending") ///
 mlcolor(black ..) msize(*1.5 ..) ///
 graphregion(color(white)) bgcolor(white) ///
 legend(pos(11) ring(0) col(1))

graph export `"${PATH_OUT_FIGURES}/ms_adv.pdf"', replace

* use the Normalized ad spending

scatter mshare_yrl? ad_spending_norm, ms(O ..) ///
 ytitle("Average yearly market share") ///
 xtitle("Yearly advertisement spending by customer") ///
 mlcolor(black ..) msize(*1.5 ..) ///
 graphregion(color(white)) bgcolor(white) ///
 legend(ring(0) bplacement(neast)  col(1))

graph export `"${PATH_OUT_FIGURES}/ms_adv_quick.pdf"', replace


restore

*******************************************************************************
* Awareness indicator
preserve

*drop trash
keep month vtest

*drop duplicates
qui by month vtest, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*prep data
sort month
gen newdate=date(string(month),"YM")
gen date=newdate
format %td date	

twoway line vtest date, lcolor(cranberry) ///
	xline(18993 19359 19724 20089 20454, lcolor(bl) lpattern(dash) lwidth(thin)) ///
	xtitle("") ///
	graphregion(color(white)) bgcolor(white) ///
	ytitle("Share of fully-informed")
	
graph export `"${PATH_OUT_FIGURES}/awareness.pdf"', replace

restore

*******************************
* Awareness and switching rate - graph and correlation
preserve

*drop trash
keep month vtest switching

*drop duplicates
qui by month vtest, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*prep data
sort month
gen newdate=date(string(month),"YM")
gen date=newdate
format %td date	

label var vtest "Share of fully-informed"
twoway (line vtest date, yaxis(1)) (line switching date, lwidth(thick) yaxis(2)), ///
	xline(18993 19359 19724 20089 20454, lcolor(bl) lpattern(dash) lwidth(thin)) ///
	ytitle("Share of fully informed consumers",axis(1)) ///
	ytitle("Churn rate",axis(2)) ///
	xtitle("") ///	
	graphregion(color(white)) bgcolor(white) ///
	legend(label(1 Information indicator) label(2 Supplier churn rate))	
		
	
graph export `"${PATH_OUT_FIGURES}/awareness_switching.pdf"', replace
* Export data for redoing graph in python.
export delimited using `"${PATH_OUT_DATA}/awareness_switching_data.csv"', replace 
restore

* is there a correlation?
preserve
keep month vtest switching 

twoway (scatter vtest switching)
gen vtest_minus1=vtest[_n-1]
gen vtest_minus2=vtest[_n-2]
gen vtest_minus3=vtest[_n-3]
twoway (scatter vtest_minus3 switching)

gen year=int(month/100)
*drop duplicate
qui by year switching, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup month

*Avg vtest
by year, sort: egen avg_vtest=mean(vtest)
*total switch
by year, sort: egen total_switch=sum(switching)

drop if year==2016 //no complete info (only until half of the year

twoway (scatter avg_vtest total_switch)

restore

*******************************
* Awareness and advertisement
preserve

*drop trash
keep month vtest switching firm ad_spending_mly

*create monthly figures (no difference by supplier)
by month, sort: egen ad_spending=total(ad_spending_mly)

*drop duplicates
qui by month ad_spending, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*prep data
sort month
gen newdate=date(string(month),"YM")
gen date=newdate
format %td date	

label var vtest "Share of fully-informed"
label var ad_spending "Advertisement"

twoway (line vtest date, lwidth(thick)) (line ad_spending date, yaxis(2)), ///
	xline(18993 19359 19724 20089 20454, lcolor(bl) lpattern(dash) lwidth(thin)) ///
	xtitle("") ///
	ytitle("") ///	
	graphregion(color(white)) bgcolor(white) ///
	legend(label(1 Fully-informed) label(2 Ad spending (in 1,000 EUR)) )	
		
graph export `"${PATH_OUT_FIGURES}/awareness_advertisement.pdf"', replace

restore

*******************************************************************************
* PRICES by month
preserve

*drop trash
keep firm month contract_type price mshare mshare_contract

*prep data
sort month
gen newdate=date(string(month),"YM")
gen date=newdate
format %td date	

label var price "Price"
label var mshare_contract "Market Share"

*graph
*ONLY CONV
twoway (line price date if firm=="ECS" & contract_type==1, color(navy)) ///
	(line price date if firm=="EDF" & contract_type==1, color(maroon)) ///
	(line price date if firm=="Eni" & contract_type==1, color(dkorange)) ///
	(line price date if firm=="Essent" & contract_type==1, color(gray)) ///
	(line price date if firm=="Other", color(sand)), ///
	yscale(range(20 40)) ///
	xtitle("") ///
	ytitle("Average real price")	///
	title("Conventional contracts") ///
	legend(label(1 ECS) label(2 EDF) label(3 Eni) label(4 Essent) label(5 Outside)) ///
	graphregion(color(white)) bgcolor(white) ///
	name(conventional,replace)  

graph export `"${PATH_OUT_FIGURES}/conv_prices.pdf"', replace
* Export data for redoing graph in python.
export delimited using `"${PATH_OUT_DATA}/prices_data.csv"', replace 
*ONLY GREEN
twoway (line price date if firm=="ECS" & contract_type==3, color(navy)) ///
	(line price date if firm=="EDF" & contract_type==3, color(maroon)) ///
	(line price date if firm=="Eneco", color(green)) ///
	(line price date if firm=="Eni" & contract_type==3, color(dkorange)) ///
	(line price date if firm=="Essent" & contract_type==3, color(gray)) ///
	(line price date if firm=="Lampiris", color(cranberry)), ///
	yscale(range(20 40)) ///
	xtitle("") 	///
	ytitle("Average real price")	///
	title("Green contracts") ///
	legend(label(1 ECS) label(2 EDF) label(3 Eneco) label(4 Eni) label(5 Essent) label(6 Lampiris)) ///
	graphregion(color(white)) bgcolor(white) ///
	name(green,replace) 
graph export `"${PATH_OUT_FIGURES}/green_prices.pdf"', replace