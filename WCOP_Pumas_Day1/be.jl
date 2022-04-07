
# Call Other Necessary Packages 
using CSV
using DataFramesMeta
using Pumas
using PumasUtilities
using Bioequivalence 
#using Bioequivalence.GLM: lm, @formula





# Load Data
## This is a 2x2 nonreplicated crossover design (RT | TR) with 2 periods
## Dataframe must have subject ID, their respective sequence, period, and exposure of interest (AUC Cmax) 
df_be = CSV.read("WCOP//be_data//be_data.csv", DataFrame, missingstrings=["NA", ".", ""])




# Run Bioequivalence analysis 
## endpoint - indicate exposure measure you want to use to compare for BE 
## method   - indicate your statistical method 
##           :fda (default value; choose linear models if nonreplicated and linear mixed models if replicated)
##           :lm use a linear model
##           :lmm use a linear mixed model
##           :nonpar use a nonparametric model (default if name of endpoint includes tmax ignoring case)
## reml - choosing between restricted maximum likelihood objective (reml) and maximum likehihood objective (ml)
##      - select reml if we want to match SAS default values 

output_cmax = pumas_be(df_be, endpoint = :Cmax, method = :fda, reml = true)
output_auc  = pumas_be(df_be, endpoint = :AUC, method = :fda, reml = true)



# Interpreting results 
## output in REPL will state the design, sequence, periods, and subjects per sequence 
## output table:
##       - PE = point estimate for the natural log of the geometric mean ratio 
##       - SE = standard error 
##       - lnLB = natural log of the lower bound of the geometric mean ratio 
##       - lnUB = natural log of the upper bound of the geometric mean ratio 
##       - GMR = geometric mean ratio 
##       - LB = lower bound of the GMR 
##       - UB = upper bound of the GMR 

## can use these results to compare to bioequivalence criteria 
## If we pick our criteria to be 80% to 125%, both Cmax and AUC pass 
## If we pick our criteria to be between 90% and 110%, only AUC passes 



# How to look at individual outputs:
## comparing formulations (R vs T) 
output_auc.data_stats.formulation

## comparing sequences (RT vs TR)
output_auc.data_stats.sequence

## comparing periods (1 vs 2)
output_auc.data_stats.period

## output statistical model results 
output_auc.model

## perform Wald test - assesses whether the model parameters are jointly statistically significant from zero
output_auc.model_stats.Wald

## compare least squares geometric means for the formulation 
output_auc.model_stats.lsmeans

