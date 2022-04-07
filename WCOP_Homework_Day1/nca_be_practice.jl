
# Call Necessary Packages 
using CSV
using Chain
using DataFrames
using Dates
using NCA
using NCAUtilities
using NCA.Unitful
using PumasUtilities




###################################################################
#                                                                 #
#        PART 1 - NON-COMPARTMENTAL ANALYSIS PRACTICE             #
#                                                                 #
###################################################################

# Read in nca_homework.csv file under the "data" folder of WCOP Homework as a dataframe
## Look at your data in the REPL or the Julia environment 
## What route of administration is indicated? Is this a case of multiple or single dose? 


# Map your dataframe to a NCA Population
## Are there multiple doses? If so, what do you need to map to ensure you can compare between doses?



# Generate a dataframe of individual NCA Parameters for AUC from 0 to 24, Cmax (not dose-normalized), Terminal volume, Clearance, and Half-life



# Run an NCA that includes report formatting data 
## Study ID: STUDY-HOMEWORK
## Study title: "NCA Homework Practice" 
## Author: Your name 
## Sponsor: WCOP Pumas Workshop 



# Indicate what parameters you would like to summarize in your final report 
## The parameters to include should be tmax, cmax, aucinf_obs, vz_f_obs, and cl_f_obs


# Generate an NCA report and save the PDF 



# Does this drug exhibit dose proportionality? 











###################################################################
#                                                                 #
#              PART 2 - BIOEQUIVALENCE ANALYSIS                   #
#                                                                 #
###################################################################


# Read in be_homework.csv file under the "data" folder of WCOP Homework as a dataframe




# Perform a bioequivalence analysis for both Cmax and AUC using whichever method is best for this trial design 
## Hint: Is this design replicated or nonreplicated? Based on this, which statistical method should you use? 



# Interpret your results and answer the following questions: 
## 1. What is the geometric mean ratio and associated lower/upper bounds for the BE analysis of AUC? 
## 2. What is the geometric mean ratio and associated lower/upper bounds for the BE analysis of Cmax? 
## 3. If our criteria for bioequivelance is to acheive a GMR between 80% and 125%... 
##            a. Does the tablet & capsule form acheive bioequivalence based on AUC? 
##            b. Does the tablet & capsule form acheive bioequivalence based on Cmax? 
## 4. If our criteria for bioequivelance is to acheive a GMR between 90% and 110%... 
##            a. Does the tablet & capsule form acheive bioequivalence based on AUC? 
##            b. Does the tablet & capsule form acheive bioequivalence based on Cmax? 

