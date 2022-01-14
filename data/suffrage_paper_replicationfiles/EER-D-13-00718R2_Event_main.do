clear

*The next line reads in the data for the baseline reform measure. Note that the path should be corrected.

*use "........\EER-D-13-00718R2_Event_study_main_sample_reform.dta"

*Table 8 except for regression 7 which estimated below 

xi: logit Reform1  Revolution2_2012 laggdp   lagn  lagurban  war ww1   lagenrollment Gold_standard timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  Revolution_ling_2012 laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  lagrevolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b New_learn_ling laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: relogit Reform1  revolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: xtlogit Reform1  revolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, fe

* Table 10 event study part

xi: logit Reform1  revolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard lagtrade  timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard lagtc  timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard lagagricultural  timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b laggdp   lagn  lagurban  warintensity2 ww1  lagenrollment Gold_standard timesincereform _spline*, cluster(ccode) 


* Table D1 in the supplementary appendix


xi: logit Reform1  revolution5b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b trend_hp cycle_hp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b g0_rev laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

xi: logit Reform1  revolution4b g2_3rev laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 

*Table 8, regression 7
 
clear

*The next line reads in the data for the reform measure which allows for reversals. Note that the path should be corrected.

*use "........\EER-D-13-00718R2_Event_study_main_sample_reversals.dta"


xi: logit Reformreversal revolution4b laggdp   lagn  lagurban  war ww1  lagenrollment Gold_standard   timesincereform _spline*, cluster(ccode) 
