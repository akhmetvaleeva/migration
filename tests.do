set more off
import delimited "di.csv", encoding(UTF-8) clear 
xtset id year

xtreg lcoefmigr lpop ipc davto dlong destprir infmort dlincome goods gini belowmin incgroup dlbadwater dlvrppc percmin dtotsp famsubs docload unemp

/*
foreach v of varlist lpribyv lvybyv lcoefmigr lpop ipc avto v13 estprir infmort lincome goods gini belowmin incgroup lbadwater lvrppc percmin totsp famsubs docload unemp { 
    xtunitroot ht `v'
	xtunitroot llc `v'
	xtunitroot ips `v'
}
*/

spatwmat using "dij3.dta", name(W1) standardize
spatwmat using "dij4.dta", name(W2) standardize
spatwmat using "neighb2.dta", name(W3) binary standardize

set more off
forvalues i = 2006(1)2014 {
preserve
drop if year != `i'
spatgsa lcoefmigr, weights(W1) moran
spatgsa lcoefmigr, weights(W2) moran
spatgsa lcoefmigr, weights(W3) moran
restore
}

* net install xsmle, all from(http://www.econometrics.it/stata)
set more off
local xs " lpop ipc davto dlong destprir infmort dlincome goods gini belowmin incgroup dlbadwater dlvrppc percmin dtotsp famsubs docload unemp"
* BFGS converges better
xsmle lcoefmigr `xs', wmat(W1) fe type(both) model(sar) cluster(id) hausman
mat beta=e(b)
est store sar1
xsmle lcoefmigr `xs', wmat(W2) fe type(both) model(sar) cluster(id)   hausman
est store sar2
xsmle lcoefmigr `xs', wmat(W3) fe type(both) model(sar) cluster(id) hausman
est store sar3
xsmle lcoefmigr `xs', emat(W1) fe type(both) model(sem) cluster(id) hausman
est store sem1
xsmle lcoefmigr `xs', emat(W2) fe type(both) model(sem) cluster(id) hausman
est store sem2
xsmle lcoefmigr `xs', emat(W3) fe type(both) model(sem) cluster(id) hausman
est store sem3

* ssc install outreg2
outreg2 [sar1 sar2 sar3] using spatial.xls, excel replace
outreg2 [sem1 sem2 sem3] using spatial2.xls, replace

* Without Moscow and SpB
spatwmat using "dij_nomoscow.dta", name(Wm) standardize
preserve
drop if id==78 | id==79
local xs " lpop ipc davto dlong destprir infmort dlincome goods gini belowmin incgroup dlbadwater dlvrppc percmin dtotsp famsubs docload unemp"
xsmle lcoefmigr `xs' , wmat(Wm) fe type(both) model(sar) cluster(id) hausman
restore

