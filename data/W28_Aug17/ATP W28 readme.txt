PEW RESEARCH CENTER
Wave 28 American Trends Panel 
Dates: August 8-21, 2017
Mode: Web
Language: English and Spanish
N=4,971

***************************************************************************************************************************
NOTES

The reports “Wide Partisan Gaps in U.S. Over How Far the Country Has Come on Gender Equality” and “On Gender Differences, No Consensus on Nature vs. Nurture” and some related 
Fact Tank posts used questions from both Wave 28 and Wave 29 and used only the respondents who responded to both waves. In order to replicate the analysis found in these 
publications, the weight WEIGHT_W28W29 must be applied and the data must be filtered on the variable W28W29_FLAG=1 to identify respondents to both Wave 28 and Wave 29. 
The datasets for W28 and W29 may also be combined.


There was a programming error for question SOCDISCW – the word ‘believe’ was missing. A total of 1,085 cases were shown the incorrect version. A correct version 
was re-asked as a standalone question to all N=1,085 respondents who got the incorrect wording. A total of N=937 responded, the remaining N=148 were coded in the 
data as SKIPPED. The affected cases are flagged in the data in variable SOCDISCW_FLAG_W28.


In order to filter data to show only never-married Americans like in the September 14, 2017 blog post “As U.S. marriage rate hovers at 50%, 
education gap in marital status widens” (http://www.pewresearch.org/fact-tank/2017/09/14/as-u-s-marriage-rate-hovers-at-50-education-gap-in-marital-status-widens/), 
use the filter F_MARITAL_FINAL=6 OR (F_MARITAL_FINAL=2 AND LWP_RF2=2).


The following variables are included in the dataset as coded open end response. Responses to the first three mentions have been included. 
SOCVALM_m1_W28
SOCVALM_m2_W28
SOCVALM_m3_W28
SOCDISCM_m1_W28
SOCDISCM_m2_W28
SOCDISCM_m3_W28
SOCVALW_m1_W28
SOCVALW_m2_W28
SOCVALW_m3_W28
SOCDISCW_m1_W28
SOCDISCW_m2_W28
SOCDISCW_m3_W28

***************************************************************************************************************************
WEIGHTS 

The W28 dataset contains three weight variables:

WEIGHT_W28 is the weight for the standalone W28 respondents except for questions NEWS_PLATFORM, SNS, and SNSNEWS.

WEIGHT_W28_SOCIALMEDIA is the weight for the standalone W28 respondents but is used only for NEWS_PLATFORM, SNS, and SNSNEWS. This weight rakes to news users from social media site. 
This is the same process that was used in W14. Topline numbers and numbers used in reports for these questions use this weight. See W28 methodology for more details on raking targets.

WEIGHT_W28W29 is a longitudinal weight for N=4,573 respondents to both W28 and W29.





