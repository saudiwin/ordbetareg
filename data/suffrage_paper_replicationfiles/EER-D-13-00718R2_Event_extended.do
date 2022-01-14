clear

*The next line reads in the data for the baseline reform indicator. Note that the path should be corrected.

*use "......\EER-D-13-00718R2_Event_study_extended_sample_reform.dta"

xi: logit reform1 Revolution4b war ww1 timesincereformbase _splinebase*, cluster(ccode)

xi: relogit reform1 Revolution4b war ww1 timesincereformbase _splinebase*, cluster(ccode)

xi: xtlogit reform1 Revolution4b war ww1 timesincereformbase _splinebase*, fe 

xi: logit reform1 Revolution4b laggdp lagn war ww1 timesincereformbase _splinebase*, cluster(ccode)

xi: relogit reform1 Revolution4b laggdp lagn war ww1 timesincereformbase _splinebase*, cluster(ccode)

xi: xtlogit reform1 Revolution4b laggdp lagn war ww1 timesincereformbase _splinebase*, fe

*Equation (16) -- Linear probability models with year effects

xi: reg reform1 Revolution4b war i.ccode i.year,  cluster(ccode)

xi: reg reform1 Revolution4b laggdp lagn war i.ccode i.year,  cluster(ccode)


* Models with reversals (columns 4 and 8 in Table 9)

clear

*The next line reads in the data for the reform indicator which allows for reversals. Note that the path should be corrected.

*use "......\EER-D-13-00718R2_Event_study_extended_sample_reversals.dta"

xi: logit reform_reversal Revolution4b war ww1 timesincereformreversal _splinereversal*, cluster(ccode)

xi: logit reform_reversal Revolution4b laggdp lagn war ww1 timesincereformreversal _splinereversal*, cluster(ccode)



