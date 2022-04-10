
## SUMMARY OF TODAY'S POPPK WORKFLOW 
# 1. Read in the givwn CSV file into a DataFrame 
# 2. Gain an understanding of your data & initial estimates by performing a simple NCA 
# 3. Base model building (residual error model, # of compartments, etc.)
# 4. Perform covariate analysis & include appropriate covariates into model 
# 5. Simulation with final model selected 


using CSV
using DataFramesMeta
using Dates
using Pumas
using Pumas.Latexify
using PumasUtilities
using Random
using CairoMakie
using Serialization


# Read in data
## iv infusion given over 2 hours with demographic information (age, weight, sex, crcl)

pkdata = CSV.read("poppk_data/nlme_sample.csv", DataFrame; missingstrings=["NA",""])


## filter only first dose data for NCA preview 
pkdata_nca = @rsubset pkdata :OCC == 1


## read in nca data 
pop_nca = read_nca(pkdata_nca,
                id            = :ID,
                time          = :TIME,
                observations  = :DV,
                amt           = :AMT,
                route         = :ROUTE,
                duration      = :DURATION,
                group         = [:GROUP,],
                llq           = 0.001)
#

# Preview Data 
## individual plot - linear scale 
obsvstimes = observations_vs_time(pop_nca[1])

## individual plot - semi-log scale
obsvstimes = observations_vs_time(pop_nca[1], axis = (yscale = log,))

## Grid of individual plots
ctplots = observations_vs_time(pop_nca[1:9], 
                                axis = (xlabel = "Time (hour)", 
                                        ylabel = "Drug Concentration (mg/L)",
                                        yscale = log,),
                                paginate = true, #creates multiple pages  
                                columns = 3, rows = 3, #number of col/rows per page 
                                facet = (combinelabels = true,)) #creates 1 label for each page
ctplots[1] # page 1 

## Mean concentration-time curve of population 
summary_observations_vs_time(pop_nca,
                                axis = (xlabel = "Time (hour)", 
                                ylabel = "Drug Concentration (mg/L)"))
#                        

## Run simple NCA 
nca_for_poppk = run_nca(pop_nca, sigdigits=3)

## Summarize for initial estimates / initial understanding of PK 
param_summary  = summarize(nca_for_poppk.reportdf, 
                                stratify_by=[:GROUP,], # stratifying by group so we can compare each dose 
                                parameters = [:vz_obs, 
                                          :cl_obs, 
                                          :aucinf_obs,
                                          :cmax])
#





##########################################################
#                                                        #
#                   1 COMPRATMENT MODEL                  #
#                                                        #
##########################################################


pop = read_pumas(pkdata,
                  id=:ID,
                  time=:TIME,
                  amt=:AMT,
                  covariates = [:WT,:AGE,:SEX,:CRCL,:GROUP],
                  observations = [:DV],
                  cmt=:CMT,
                  evid=:EVID,
                  rate=:RATE)
#



mdl_1cmp_comb = @model begin
    @metadata begin
      desc = "base model: 1comp"
      timeu = u"hr"
    end 

    @param begin
      "Clearance (L/hr)"
      tvcl   ∈ RealDomain(lower = 0.0001)
      "Volume (L)"
      tvvc   ∈ RealDomain(lower = 0.0001)
      """
      - ΩCL
      - ΩVc
      """
      Ω      ∈ PDiagDomain(2) #covariance matrix 
      # σ_add  ∈ RealDomain(lower = 0.0001) # standard deviation 
      # σ_prop ∈ RealDomain(lower = 0.0001) # standard deviation 
      "Additive RUV"
      σ²_add ∈ RealDomain(lower = 0.0001) # variance 
      "Proportional RUV"
      σ²_prop ∈ RealDomain(lower = 0.0001) # variance  
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @covariates WT AGE SEX CRCL GROUP 

    @pre begin
      CL = tvcl * exp(η[1])
      Vc = tvvc * exp(η[2])
    end

    @dynamics Central1

    @derived begin
      CONC := @. Central/Vc #suppressing CONC from being output using ":="
      # additive error model 
      #DV ~ @. Normal(CONC, σ_add) # using standard deviation 
      #DV ~ @. Normal(CONC, sqrt(σ²_add)) # using variance 
      # proportional error model 
      #DV ~ @. Normal(CONC, abs(CONC)*σ_prop) # using standard deviation 
      #DV ~ @. Normal(CONC, sqrt(CONC^2*σ²_prop))   # using variance 
      # combination error model 
      """
      DrugY Concentration (ng/mL)
      """
      DV ~ @. Normal(CONC, √(abs(CONC)^2*σ²_prop + σ²_add)) # using variance
    end
end



params_1cmp_comb = (
    tvvc = 5,
    tvcl = 0.2,
    Ω = Diagonal([0.01, 0.01]),
    σ²_add = 0.01,
    σ²_prop = 0.01)
#

## Initial estimate explorer
ee_1cmp_comb = explore_estimates(mdl_1cmp_comb, 
                                 pop, 
                                 params_1cmp_comb)
# generate latex code in the REPL
latexify(mdl_1cmp_comb, :dynamics)
# render the latex code in the plot pane
# you can right click on the rendered code and save the tex to use in Word e.g.
render(latexify(mdl_1cmp_comb, :dynamics))

## Update coefficients based on app
params_1cmp_comb = coef(ee_1cmp_comb)

## Create a DataFrame from the exploration
DataFrame(ee_1cmp_comb)

## Evaluate initial loglikelihood 
pkfit_1cmp_comb = loglikelihood(mdl_1cmp_comb,
                                pop,
                                params_1cmp_comb,
                                Pumas.FOCE())

## explore influential individuals 
pkfit_1cmp_comb = findinfluential(mdl_1cmp_comb,
                                  pop,
                                  params_1cmp_comb,
                                  Pumas.FOCE())                
## Run model fitting 
pkfit_1cmp_comb = fit(mdl_1cmp_comb,
                      pop,
                      params_1cmp_comb,
                      Pumas.FOCE())
#


## Save your fitted model 
serialize("poppk_data/base_model_building/pkfit_1cmp_comb.jls", pkfit_1cmp_comb)

## Call your saved fitted model 
pkfit_1cmp_comb = deserialize("poppk_data/base_model_building/pkfit_1cmp_comb.jls")


## Model diagnostics to compare models - but there is also another way using evaluate_diagnostics (shown next)
metrics_pkfit_1cmp_comb = metrics_table(pkfit_1cmp_comb) 

# prtint the coefficients of the model fitting
coefficients_table(pkfit_1cmp_comb)

## Run some post-processing steps that can be used later to compare our different compartment models 
infer_1cmp_comb = infer(pkfit_1cmp_comb)
coeftable(infer_1cmp_comb) # Get infer output as a dataframe
infer_1cmp_comb.fpm # The fitted model summary used for inference
infer_1cmp_comb.level # confidence level
infer_1cmp_comb.vcov # variance-covariance matrix

inspect_1cmp_comb = inspect(pkfit_1cmp_comb)
goodness_of_fit(inspect_1cmp_comb)
dfinsp = DataFrame(inspect_1cmp_comb)
DataFrame(inspect_1cmp_comb, include_events = false)
DataFrame(inspect_1cmp_comb, include_covariates = false)

inspect_1cmp_comb_FO = inspect(pkfit_1cmp_comb, wres_approx = Pumas.FO())
DataFrame(inspect_1cmp_comb_FO, include_covariates = false, include_events = false)

inspect_1cmp_comb_npde = inspect(pkfit_1cmp_comb, nsim = 100, rng = Random.seed!(1234))
DataFrame(inspect_1cmp_comb_npde, include_covariates = false, include_events = false)

##########################################################
#                                                        #
#                   2 COMPRATMENT MODEL                  #
#                                                        #
##########################################################

mdl_2cmp_comb = @model begin
    @param begin
      tvcl   ∈ RealDomain(lower = 0.0001)
      tvvc   ∈ RealDomain(lower = 0.0001)
      tvq    ∈ RealDomain(lower = 0.0001)
      tvvp   ∈ RealDomain(lower = 0.0001)
      Ω      ∈ PDiagDomain(2)
      σ²_add  ∈ RealDomain(lower = 0.0001)
      σ²_prop ∈ RealDomain(lower = 0.0001)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @covariates WT AGE SEX CRCL GROUP

    @pre begin
      CL = tvcl * exp(η[1])
      Vc = tvvc * exp(η[2])
      Q  = tvq  
      Vp = tvvp 
    end

    @dynamics Central1Periph1

    @derived begin
      CONC := @. Central/Vc
      DV ~ @. Normal(CONC, √(abs(CONC)^2*σ²_prop + σ²_add)) # using variance
    end
end



params_2cmp_comb = (
    tvvc   = 5,
    tvcl   = 0.02,
    tvq    = 0.01,
    tvvp   = 10,
    Ω      = Diagonal([0.01, 0.01]),
    σ²_add  = 0.01,
    σ²_prop = 0.01)
#



## Maximum likelihood estimation
pkfit_2cmp_comb = fit(mdl_2cmp_comb,
                      pop,
                      params_2cmp_comb,
                      Pumas.FOCE())
#


# Serialize fits 
serialize("poppk_data/base_model_building/pkfit_2cmp_comb.jls", pkfit_2cmp_comb)
pkfit_2cmp_comb = deserialize("poppk_data/base_model_building/pkfit_2cmp_comb.jls")


## Model diagnstic - can be used to compare different compartment models manually 
metrics_pkfit_2cmp_comb = metrics_table(pkfit_2cmp_comb) 

## Compare the estimates with 1-cmp model
compare_estimates(;pkfit_2cmp_comb, pkfit_1cmp_comb)

## Create dataframe to compare for base model selection 
df_compartment_comp = innerjoin(metrics_pkfit_2cmp_comb, metrics_pkfit_1cmp_comb, on=:Metric, makeunique=true)

#


# Run some post-processing steps to compare 1 cmt and 2 cmt 
infer_2cmp_comb = infer(pkfit_2cmp_comb)
inspect_2cmp_comb = inspect(pkfit_2cmp_comb) 



## can also use to look at Compare 1-comp and 2-comp models using GUI! 
compare_basemodels = evaluate_diagnostics([
                                            (pkfit_1cmp_comb, 
                                            infer_1cmp_comb, 
                                            inspect_1cmp_comb), 
                                            
                                            (pkfit_2cmp_comb, 
                                            infer_2cmp_comb, 
                                            inspect_2cmp_comb)],

                                            categorical =[:SEX, :GROUP])
##










##########################################################
#                                                        #
#             BASE MODEL POST-PROCESSING                 #
#                                                        #
##########################################################


# selecting 2 compartment model as base model (saving fitted model as such)
serialize("poppk_data/base/pkfit_base.jls", pkfit_2cmp_comb)
pkfit_base = deserialize("poppk_data/base/pkfit_base.jls")

## post-hoc concentration predictions 
pred_base = predict(pkfit_base, obstimes=0:0.5:168)
pred_base_df = DataFrame(pred_base)
DataFrame(pred_base, include_covariates = false, include_events = false)

CSV.write("poppk_data/base/pred_base.csv", pred_base)

## individual predictions
indiv_base = reduce(vcat, DataFrame.(icoef(pkfit_base)))
CSV.write("poppk_data/base/indiv_base.csv", indiv_base)

## empirical bayes 
ebes_base = DataFrame(empirical_bayes(pkfit_base))
empirical_bayes_dist(inspect_2cmp_comb)
CSV.write("poppk_data/base/ebes_base.csv", ebes_base)

## residuals 
resid_base = DataFrame(wresiduals(pkfit_base))
CSV.write("poppk_data/base/resid_base.csv", resid_base)

## individual subject fits
sf = subject_fits(inspect_2cmp_comb, 
              separate = true, 
              columns = 2, rows = 2, paginate = true,
              facet = (combinelabels = true,),
              figure = (resolution = (1400, 1000),fontsize = 36),
              axis = (ylabel = "Observed\nPredicted drugY (ng/mL)",
              yscale=log10, ytickformat=x -> string.(round.(x; digits=1)), 
                                            ygridwidth = 3, 
                                            yminorgridcolor = :darkgrey,
                                            yminorticksvisible = true,
                                            yminorgridvisible = true,
                                            yminorticks = IntervalsBetween(10),
                                            xminorticksvisible = true,
                                            xminorgridvisible = true,
                                            xminorticks = IntervalsBetween(5),
                                            # limits = (nothing, nothing, nothing, 30),
                                            spinewidth = 2))  
sf[2]

#########################################################################
#                                                                       #
#                      COVARIATE MODEL BUILDING                         #
#                                                                       #
#########################################################################

pop = read_pumas(pkdata,
                    id=:ID,
                    time=:TIME,
                    amt=:AMT,
                    covariates=[:WT,:AGE,:SEX,:CRCL,:GROUP],
                    observations=[:DV],
                    cmt=:CMT,
                    evid=:EVID,
                    rate=:RATE)
#

## some basic plots
covariates_check(pop) ### static plot
covariates_check(pop, covariates =[:WT]) |> interactive #interactive plot
covariates_dist(pop)
observations_vs_time(pop)
observations_vs_time(pop, axis = (yscale = log10,))


###################### WEIGHT ##########################################


mdl_base_wt = @model begin
    @param begin
      tvcl   ∈ RealDomain(lower = 0.0001)
      tvvc   ∈ RealDomain(lower = 0.0001)
      tvq    ∈ RealDomain(lower = 0.0001)
      tvvp   ∈ RealDomain(lower = 0.0001)
      Ω      ∈ PDiagDomain(2)
      σ²_add  ∈ RealDomain(lower = 0.0001)
      σ²_prop ∈ RealDomain(lower = 0.0001)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @covariates WT 

    @pre begin
        wtCL = (WT/70)^0.75
        wtV  = (WT/70)
        CL   = tvcl * wtCL * exp(η[1])
        Vc   = tvvc * wtV  * exp(η[2])
        Q    = tvq  
        Vp   = tvvp 
    end

    @dynamics Central1Periph1

    @derived begin
      CONC := @. Central/Vc
      DV ~ @. Normal(CONC, √(abs(CONC)^2*σ²_prop + σ²_add)) # using variance
    end
end



params_base_wt = (
      tvvc   = 5,
      tvcl   = 0.02,
      tvq    = 0.01,
      tvvp   = 10,
      Ω      = Diagonal([0.01, 0.01]),
      σ²_add  = 0.01,
      σ²_prop = 0.01)
#


## Maximum likelihood estimation
pkfit_base_wt = fit(mdl_base_wt,
              pop,
              params_base_wt,
              Pumas.FOCE())
#

serialize("poppk_data/base_wt/pkfit_base_wt.jls", pkfit_base_wt)
#pkfit_base_wt = deserialize("./templates/data/base_wt/pkfit_base_wt.jls")


## model diagnstics
metrics_pkfit_base_wt = metrics_table(pkfit_base_wt) 
df_comp_base_wt = innerjoin(metrics_pkfit_base_wt, metrics_pkfit_2cmp_comb, on=:Metric, makeunique=true)


## post-processing 
inspect_base_wt = inspect(pkfit_base_wt)
empirical_bayes_vs_covariates(inspect_base_wt)
goodness_of_fit(inspect_base_wt)

###################### CRCL ##########################################


mdl_base_wt_crcl = @model begin
    @param begin
      tvcl   ∈ RealDomain(lower = 0.0001)
      tvvc   ∈ RealDomain(lower = 0.0001)
      tvq    ∈ RealDomain(lower = 0.0001)
      tvvp   ∈ RealDomain(lower = 0.0001)
      #tvβ    ∈ RealDomain(lower = 0.0001)
      Ω      ∈ PDiagDomain(2)
      σ²_add  ∈ RealDomain(lower = 0.0001)
      σ²_prop ∈ RealDomain(lower = 0.0001)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @covariates WT CRCL 

    @pre begin
        wtCL     = (WT/70)^0.75
        wtV      = (WT/70)
        crcl_eff = (CRCL/95)^0.75
        #β        = tvβ #* exp(η[3])
        #crcl_eff = (CRCL/95)^β
        CL       = tvcl * crcl_eff * wtCL * exp(η[1])
        Vc       = tvvc * wtV  * exp(η[2])
        Q        = tvq  
        Vp       = tvvp 
    end

    @dynamics Central1Periph1

    @derived begin
      CONC := @. Central/Vc
      DV ~ @. Normal(CONC, √(abs(CONC)^2*σ²_prop + σ²_add)) # using variance
    end
end



params_base_wt_crcl = (
        tvvc   = 5,
        tvcl   = 0.02,
        tvq    = 0.01,
        tvvp   = 10,
        #tvβ    = 0.75,
        Ω      = Diagonal([0.01, 0.01]),
        σ²_add  = 0.01,
        σ²_prop = 0.01)
#

## Maximum likelihood estimation
pkfit_base_wt_crcl = fit(mdl_base_wt_crcl,
              pop,
              params_base_wt_crcl,
              Pumas.FOCE())
#

serialize("poppk_data/base_wt_crcl/pkfit_base_wt_crcl.jls", pkfit_base_wt_crcl)
#pkfit_base_wt_crcl = deserialize("./templates/data/base_wt_crcl/pkfit_base_wt_crcl.jls")

## list all models that are in the environment
list_models()

## model diagnstics
metrics_pkfit_base_wt_crcl = metrics_table(pkfit_base_wt_crcl)
df_comp_base_wt_crcl = innerjoin(metrics_pkfit_base_wt_crcl, metrics_pkfit_base_wt, metrics_pkfit_2cmp_comb, on=:Metric, makeunique=true)


## post-processing 
inspect_base_wt_crcl = inspect(pkfit_base_wt_crcl)
empirical_bayes_vs_covariates(inspect_base_wt_crcl)
goodness_of_fit(inspect_base_wt_crcl)

## Save coefficients for model parameter estimates to a text document 
open("poppk_data/base_wt_crcl/coef_base_wt_crcl.txt", "w") do io
  println(io, coeftable(pkfit_base_wt_crcl))
end


## Make an inference
infer_base_wt_crcl = infer(pkfit_base_wt_crcl)
infer_base_wt_crcl_df = coeftable(infer_base_wt_crcl)
CSV.write("poppk_data/base_wt_crcl/infer_base_wt_crcl.csv", infer_base_wt_crcl_df)

## Generate diagnostics
inspect_base_wt_crcl = inspect(pkfit_base_wt_crcl)
inspect_base_wt_crcl_df = DataFrame(inspect_base_wt_crcl)
CSV.write("poppk_data/base_wt_crcl/inspect_base_wt_crcl.csv", inspect_base_wt_crcl_df)

## post-hoc concentration predictions 
pred_base_wt_crcl = DataFrame(predict(pkfit_base_wt_crcl, obstimes=0:0.5:168))
CSV.write("poppk_data/base_wt_crcl/pred_base_wt_crcl.csv", pred_base_wt_crcl)

## individual predictions
indiv_base_wt_crcl = reduce(vcat, DataFrame.(icoef(pkfit_base_wt_crcl)))
CSV.write("poppk_data/base_wt_crcl/indiv_base_wt_crcl.csv", indiv_base_wt_crcl)

## empirical bayes 
ebes_base_wt_crcl = DataFrame(empirical_bayes(pkfit_base_wt_crcl))
CSV.write("poppk_data/base_wt_crcl/ebes_base_wt_crcl.csv", ebes_base_wt_crcl)

## residuals 
resid_base_wt_crcl = DataFrame(wresiduals(pkfit_base_wt_crcl))
CSV.write("poppk_data/base_wt_crcl/resid_base_wt_crcl.csv", resid_base_wt_crcl)

## shrinkage 
open("poppk_data/base_wt_crcl/shrinkage_base_wt_crcl.txt", "w") do io
  println(io, ηshrinkage(pkfit_base_wt_crcl))
end

## generate vpc
vpc_base_wt_crcl = vpc(pkfit_base_wt_crcl, stratify_by = [:GROUP])
vpc_plot(vpc_base_wt_crcl)


##########################################################
#                                                        #
#               EVALUATE OUR FINAL MODEL                 #
#                                                        #
##########################################################

eval_base_final = evaluate_diagnostics([(pkfit_2cmp_comb, infer_2cmp_comb, inspect_2cmp_comb), 
                                        (pkfit_base_wt_crcl, infer_base_wt_crcl, inspect_base_wt_crcl)])

## write report
report((; two_comp_model = (pkfit_base_wt_crcl, infer_base_wt_crcl, 
                            inspect_base_wt_crcl, vpc_base_wt_crcl)),
          categorical = [:SEX,:GROUP],
          date = Dates.now(),
          output = "drugy_report_2cmp", 
          clean = false, 
          title = "DrugY Population Pharmacokinetic Analysis",
          author = "Author",
          version = "v0.1",
          header = "Pumas Report",
          footer = "Confidential")


##########################################################
#                                                        #
#               FINAL MODEL SIMULATIONS                  #
#                                                        #
##########################################################


################## Simulation 1 subject ############################

# define our dosing regimen 
## giving 1000 mg at time=0 over 2 hours with an additional 5 doses given every 24 hours (total of 6 doses)
ev = DosageRegimen(1000, time=0, addl=5, ii=24, rate=500) 

# defining a subject with WT of 70 kg and CrCl of 20
s1 = Subject(id=1, events=ev, covariates=(WT=70, CRCL=20), observations=(DV=nothing,))

# defining our input parameters as the estimated parameters of our final model 
param = coef(pkfit_base_wt_crcl)

# simulating from time 1 hour to 120 hour and sampling every 1 hour (unrealistic)
subjobs = simobs(mdl_base_wt_crcl, s1, param, obstimes=1:1:120)

# can specify exact points to make simulation realistic 
# subjobs = simobs(mdl_base_wt_crcl, s1, param, obstimes=[0.001,1,2,4,8,12,24])

# convert simulated concentrations into a dataframe 
sim_df = DataFrame(subjobs) 

# plot simulated concentration 
sim_plot(mdl_base_wt_crcl,
            [subjobs], observations = [:CONC], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (mg/L)"))
#






########################### Simulate a population  ###############################

# define covaritates for WT and CRCL as being a random number between indicated ranges 
choose_covariates() = (WT = rand(55:120), 
                      CRCL = rand(20:120))
#                      
# map multiple subjects to be assigned random variables (IDs labeled as 1 to 50)
pop_for_sim = map(i->Subject(id=i, 
                            events=ev,
                            observations=(DV=nothing,),
                            covariates=choose_covariates()),
                            1:50)
#

# perform simulation 
obs_pk = simobs(mdl_base_wt_crcl, pop_for_sim, param, obstimes=0:1:120)

# convert to dataframe 
simdf = DataFrame(obs_pk)

# plotting 
sim_plot(mdl_base_wt_crcl,
            obs_pk, observations = [:CONC], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (mg/L)"))
#


sim_plot(mdl_base_wt_crcl,
            obs_pk, observations = [:DV], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (mg/L)"))
#



