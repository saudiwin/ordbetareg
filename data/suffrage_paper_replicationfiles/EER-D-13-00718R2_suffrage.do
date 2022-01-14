
clear all

*The next line reads in the data. Note that the path should be corrected.

*use "........\EER-D-13-00718R2_maindata_suffrage.dta"

*Table 2

xi: xtpcse e2c lage2c  Revolution2_2012 yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution4b yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  Revolution_ling_2012 yy* i.ccode, pairwise correlation(psar1)

*part with with controls

xi: xtpcse e2c lage2c  Revolution2_2012 laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  Revolution_ling_2012 laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode, pairwise correlation(psar1)


*Table 3

*Conley type standard errors are obtained by using x_ols available form Conley's website (see documentation file).

 gen cutoff5001=500  /*more than +-500? km away -> assumed to be independent*/
 gen cutoff5002=500  
 gen const = 1
 
xi: x_ols LatKm LongKm cutoff5001 cutoff5001 e2c lage2c revolution4b laggdp lagn lagurban lagenrollment war ww1 Gold_standard i.ccode yy*, xreg(72) coord(2)
matrix bb = e(b)
matrix v_dep = vecdiag(cov_dep)
scalar coef_e2 = bb[1,2]
scalar se_e2 = sqrt(v_dep[1,2])
display coef_e2
display coef_e2/se_e2

drop epsilon
drop dis1
drop dis2
drop window

 gen cutoff8001=800  /*more than +-800? km away -> assumed to be independent*/
 gen cutoff8002=800  

xi: x_ols LatKm LongKm cutoff8001 cutoff8002 e2c lage2c revolution4b laggdp lagn lagurban lagenrollment war ww1 Gold_standard i.ccode yy*, xreg(72) coord(2)
matrix bb = e(b)
matrix v_dep = vecdiag(cov_dep)
scalar coef_e2 = bb[1,2]
scalar se_e2 = sqrt(v_dep[1,2])
display coef_e2
display coef_e2/se_e2

drop epsilon
drop dis1
drop dis2
drop window

 gen cutoff14001=1400  /*more than +-1400? km away -> assumed to be independent*/
 gen cutoff14002=1400  
 
xi: x_ols LatKm LongKm cutoff14001 cutoff14002 e2c lage2c revolution4b laggdp lagn lagurban lagenrollment war ww1 Gold_standard i.ccode yy*, xreg(72) coord(2)
matrix bb = e(b)
matrix v_dep = vecdiag(cov_dep)
scalar coef_e2 = bb[1,2]
scalar se_e2 = sqrt(v_dep[1,2])
display coef_e2
display coef_e2/se_e2

drop epsilon
drop dis1
drop dis2
drop window
 
xi: xtpcse e2c lage2c  laggdp   lagn  lagurban  war ww1  lagrevolution4b  Gold_standard lagenrollment yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  laggdp   lagn  lagurban  war ww1  revolution4b   ownrev Gold_standard lagenrollment yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse de2  laggdp   lagn  lagurban  war ww1  lagrevolution4b Gold_standard lagenrollment yy* i.ccode, pairwise correlation(psar1)

*Table 5

xi: xtpcse e2c lage2c revolution4b ownrev laggdp   lagn  lagurban  war lagenrollment  Gold_standard  yy* i.ccode, correlation(psar1) pairwise 

xi: xtpcse e2c lage2c lagrevolution4b ownrev laggdp   lagn  lagurban  war lagenrollment  Gold_standard  yy* i.ccode, correlation(psar1) pairwise 

xi: xtpcse de2 lagrevolution4b ownrev laggdp   lagn  lagurban  war lagenrollment  Gold_standard  yy* i.ccode, correlation(psar1) pairwise 

xi: xtpcse e2c lage2c revolution4b ownrev laggdp   lagn  lagurban  war lagenrollment  Gold_standard  i.year i.ccode, het correlation(psar1) 

xi: xtpcse e2c lage2c lagrevolution4b ownrev laggdp   lagn  lagurban  war lagenrollment  Gold_standard  i.year i.ccode, het correlation(psar1) 

xi: xtpcse de2 lagrevolution4b ownrev laggdp   lagn  lagurban  war lagenrollment  Gold_standard  i.year i.ccode, het correlation(psar1) 

*Table 6 

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard g0_rev yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard g2_3rev yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution4b cycle_hp trend_hp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  year lograin laglograin revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode if year<1914, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  year revolution4b raingrowth lagraingrowth laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode if year<1914, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution5b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard yy* i.ccode if year<1914, pairwise correlation(psar1)

*regression 8 with reversals can be found at the end of this do file

*Table 7

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard yy* i.ccode, pairwise hetonly

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard yy* i.ccode, pairwise 

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard yy* i.ccode, pairwise hetonly correlation(psar1)

*Bruno correction
xi:  xtlsdvc e2c  laggdp   lagn  lagurban  war ww1  revolution4b  Gold_standard lagenrollment yy*, initial(ah) vcov(50) 

*System GMM

xi: xtdpdsys e2c laggdp   lagn  lagurban  war ww1  revolution4b   Gold_standard lagenrollment yy*, maxldep(3) maxlags(3)

*Tobit

xi: tobit e2c lage2c  laggdp   lagn  lagurban  war ww1  revolution4b Gold_standard lagenrollment yy* i.ccode, ll(0) ul(100)

*Fractional logit. Use the rescalled Suffrage variables.

xi: glm e2c_wp lage2c_wp  laggdp   lagn  lagurban  war ww1  revolution4b  Gold_standard lagenrollment yy* i.ccode if e2c_wp<=1, fam(bin) link(logit) robust cluster(ccode)

*Table 10 - suffrage part

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard lagtrade  yy* i.ccode, pairwise correlation(psar1)


xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard lagtc  yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard lagagricultural  yy* i.ccode, pairwise correlation(psar1)

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  warintensity2 ww1 lagenrollment  Gold_standard  yy* i.ccode, pairwise correlation(psar1)


*Other estimations reported or discussed in the text

*Equation (12) -- with authors of liberty

xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard authors_2 authors_3 yy* i.ccode, pairwise correlation(psar1)


*Equation (14) -- Distance interaction

xi: xtpcse e2c lage2c  laggdp   lagn  lagurban  war ww1  Revolution2_2012 Revtdist   Gold_standard lagenrollment yy* i.ccode, pairwise correlation(psar1)

*Distance interaction with year effects (discussed after Equation (14))

xi: xtpcse e2c lage2c  laggdp   lagn  lagurban  war Revtdist   Gold_standard lagenrollment i.year i.ccode, het correlation(psar1)


*Equation (15) - inequality interaction with year effects 

xi: xtpcse e2c lage2c  laggdp   lagn  lagurban  war  revolution4b lagGini Ginirev2 Gold_standard lagenrollment i.year i.ccode, hetonly correlation(psar1)

*Inequality interaction with two-year fixed effects and spatial correlation

xi: xtpcse e2c lage2c  laggdp   lagn  lagurban  war ww1  revolution4b lagGini Ginirev2 Gold_standard lagenrollment yy* i.ccode, pairwise correlation(psar1)


*Equation 3 in Supplementary Appendix D and discussed in Section 5.2.4

xi: reg e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode, cluster(year)

*Error Correction Model estimates reported in Supplementary Appendix D and discussed in Section 5.5

xi: xtpcse d.e2c lage2c  l.laggdp   l.lagn  l.lagurban  l.war l.ww1  l.revolution4b  l.Gold_standard l.lagenrollment d.laggdp   d.lagn  d.lagurban  d.war d.ww1  d.revolution4b  d.Gold_standard d.lagenrollment yy*, pairwise correlation(psar1)

clear

*use "........\EER-D-13-00718R2_reversals_suffrage.dta"

*Table 6, regression 8 which allows for democratic reversals


xi: xtpcse e2c lage2c  revolution4b laggdp   lagn  lagurban  war ww1 lagenrollment  Gold_standard  yy* i.ccode, pairwise correlation(psar1)

