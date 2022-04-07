using Pumas
using CSV
using Random
using PumasUtilities
using Serialization
using Chain 


## read in data
## iv infusion given over 2 hours with demographic information (age, weight, sex, crcl)
pkdata = DataFrame(CSV.File("WCOP//poppk_data//nlme_sample.csv", missingstrings=["NA",""]))


pop = read_pumas(pkdata,
                    id=:ID,
                    time=:TIME,
                    amt=:AMT,
                    observations=[:DV],
                    cmt=:CMT,
                    evid=:EVID,
                    rate=:RATE)
#


##########################################################
#                                                        #
#                   1 COMPRATMENT MODEL                  #
#                                                        #
##########################################################

mdl_1cmp_comb = @model begin
    @param begin
      tvcl   ∈ RealDomain(lower = 0.0001)
      tvvc   ∈ RealDomain(lower = 0.0001)
      Ω      ∈ PDiagDomain(2)
      σ_add  ∈ RealDomain(lower = 0.0001)
      σ_prop ∈ RealDomain(lower = 0.0001)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @pre begin
      CL = tvcl * exp(η[1])
      Vc = tvvc * exp(η[2])
    end

    @dynamics Central1

    @derived begin
      CONC := @. Central/Vc
      DV ~ @. Normal(CONC, (abs(CONC)*σ_prop)+σ_add)
    end
end

params_1cmp_comb = (
    tvvc = 10,
    tvcl =0.05,
    Ω = Diagonal([0.01, 0.01]),
    σ_add = 0.1,
    σ_prop = 0.1)
#

## Initial estimate explorer
ee_1cmp_comb = explore_estimates(mdl_1cmp_comb, 
                     pop, 
                     params_1cmp_comb)
#


## Update coefficients based on app
params_1cmp_comb = coef(ee_1cmp_comb)


## Run model fitting 
pkfit_1cmp_comb = fit(mdl_1cmp_comb,
                pop,
                params_1cmp_comb,
                Pumas.FOCEI(),
                ensemblealg = EnsembleThreads())
#


## Save your fitted model 
serialize("WCOP//poppk_data//pkfit_1cmp_comb.jls", pkfit_1cmp_comb)

## Call your saved fitted model 
pkfit_1cmp_comb = deserialize("WCOP//poppk_data//pkfit_1cmp_comb.jls")


## Model diagnostics to compare models 
LL_pkfit_1cmp_comb = Pumas.marginal_nll(pkfit_1cmp_comb)
minustwoLL_pkfit_1cmp_comb = 2 * LL_pkfit_1cmp_comb
AIC_pkfit_1cmp_comb = aic(pkfit_1cmp_comb)
BIC_pkfit_1cmp_comb = bic(pkfit_1cmp_comb) 













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
      Ω      ∈ PDiagDomain(4)
      σ_add  ∈ RealDomain(lower = 0.0001)
      σ_prop ∈ RealDomain(lower = 0.0001)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @pre begin
      CL = tvcl * exp(η[1])
      Vc = tvvc * exp(η[2])
      Q  = tvq  * exp(η[3])
      Vp = tvvp * exp(η[4])
    end

    @dynamics Central1Periph1

    @derived begin
      CONC := @. Central/Vc
      DV ~ @. Normal(CONC, (abs(CONC)*σ_prop)+σ_add)
    end
end



params_2cmp_comb = (
    tvvc   = 5,
    tvcl   = 0.05,
    tvq    = 0.05,
    tvvp   = 10,
    Ω      = Diagonal([0.01, 0.01, 0.01, 0.01]),
    σ_add  = 0.1,
    σ_prop = 0.1)
#



## Maximum likelihood estimation
pkfit_2cmp_comb = fit(mdl_2cmp_comb,
                pop,
                params_2cmp_comb,
                Pumas.FOCEI(),
                ensemblealg = EnsembleThreads())
#


# Serialize fits 
serialize("WCOP//poppk_data//pkfit_2cmp_comb.jls", pkfit_2cmp_comb)
pkfit_2cmp_comb = deserialize("WCOP//poppk_data//pkfit_2cmp_comb.jls")



## Model diagnstic - used to compare different compartment models  
LL_pkfit_2cmp_comb = Pumas.marginal_nll(pkfit_2cmp_comb)
minustwoLL_pkfit_2cmp_comb = 2 * LL_pkfit_2cmp_comb
AIC_pkfit_2cmp_comb = aic(pkfit_2cmp_comb)
BIC_pkfit_2cmp_comb = bic(pkfit_2cmp_comb) 




# Create dataframe to compare for base model selection 
df_compartment_comp = DataFrame(Model=["1 Compartment","2 Compartment"], 
                            AIC=[AIC_pkfit_1cmp_comb, AIC_pkfit_2cmp_comb], 
                            BIC=[BIC_pkfit_1cmp_comb, BIC_pkfit_2cmp_comb], 
                            OFV=[minustwoLL_pkfit_1cmp_comb, minustwoLL_pkfit_2cmp_comb])
#


results_eval_base = evaluate_diagnostics((; pkfit_1cmp_comb, pkfit_2cmp_comb))









##########################################################
#                                                        #
#             BASE MODEL POST-PROCESSING                 #
#                                                        #
##########################################################

# selecting 2 compartment model as base model (saving fitted model as such)
serialize("WCOP//poppk_data//base//pkfit_base.jls", pkfit_2cmp_comb)
pkfit_base = deserialize("WCOP//poppk_data//base//pkfit_base.jls")



## Save coefficients for model parameter estimates to a text document 
open("WCOP//poppk_data//base//coef_base.txt", "w") do io
  println(io, coeftable(pkfit_base))
end


## Make an inference
infer_base = infer(pkfit_base) |> coeftable
CSV.write("WCOP//poppk_data//base//infer_base.csv", infer_base)

## Generate diagnostics
inspect_base = inspect(pkfit_base) |> DataFrame
CSV.write("WCOP//poppk_data//base//inspect_base.csv", inspect_base)

## post-hoc concentration predictions 
pred_base = DataFrame(predict(pkfit_base, obstimes=0:0.5:168))
CSV.write("WCOP//poppk_data//base//pred_base.csv", pred_base)

## individual predictions
indiv_base = reduce(vcat, DataFrame.(icoef(pkfit_base)))
CSV.write("WCOP//poppk_data//base//indiv_base.csv", indiv_base)

## empirical bayes 
ebes_base = DataFrame(empirical_bayes(pkfit_base))
CSV.write("WCOP//poppk_data//base//ebes_base.csv", ebes_base)

## residuals 
resid_base = DataFrame(wresiduals(pkfit_base))
CSV.write("WCOP//poppk_data//base//resid_base.csv", resid_base)

## shrinkage 
#shrink_base = ηshrinkage(pkfit_base)
open("WCOP//poppk_data//base//shrinkage_base.txt", "w") do io
  println(io, ηshrinkage(pkfit_base))
end












#########################################################################
#                                                                       #
#                      COVARIATE MODEL BUILDING                         #
#                                                                       #
#########################################################################

pop = read_pumas(pkdata,
                    id=:ID,
                    time=:TIME,
                    amt=:AMT,
                    covariates=[:WT,:AGE,:SEX,:CRCL],
                    observations=[:DV],
                    cmt=:CMT,
                    evid=:EVID,
                    rate=:RATE)
#






###################### WEIGHT ##########################################


mdl_base_wt = @model begin
    @param begin
      tvcl   ∈ RealDomain(lower = 0.0001)
      tvvc   ∈ RealDomain(lower = 0.0001)
      tvq    ∈ RealDomain(lower = 0.0001)
      tvvp   ∈ RealDomain(lower = 0.0001)
      Ω      ∈ PDiagDomain(3)
      σ_add  ∈ RealDomain(lower = 0.0001)
      σ_prop ∈ RealDomain(lower = 0.0001)
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
        Q    = tvq  * wtCL * exp(η[3])
        Vp   = tvvp * wtV  #* exp(η[4]) removed due to overparameterization 
    end

    @dynamics Central1Periph1

    @derived begin
      CONC := @. Central/Vc
      DV ~ @. Normal(CONC, (abs(CONC)*σ_prop)+σ_add)
    end
end



params_base_wt = (
    tvvc   = 5,
    tvcl   = 0.05,
    tvq    = 0.05,
    tvvp   = 10,
    Ω      = Diagonal([0.01, 0.01, 0.01]),
    σ_add  = 0.1,
    σ_prop = 0.1)
#


## Maximum likelihood estimation
pkfit_base_wt = fit(mdl_base_wt,
              pop,
              params_base_wt,
              Pumas.FOCEI(),
              ensemblealg = EnsembleThreads())
#

serialize("WCOP//poppk_data//base_wt//pkfit_base_wt.jls", pkfit_base_wt)
#pkfit_base_wt = deserialize("./templates/data/base_wt/pkfit_base_wt.jls")


## model diagnstic 
LL_pkfit_base_wt = Pumas.marginal_nll(pkfit_base_wt)
minustwoLL_pkfit_base_wt = 2 * LL_pkfit_base_wt
AIC_pkfit_base_wt = aic(pkfit_base_wt)
BIC_pkfit_base_wt = bic(pkfit_base_wt) 


df_comp_base_wt = DataFrame(Model=["Base", "Base+WT"], 
                            AIC=[AIC_pkfit_2cmp_comb, AIC_pkfit_base_wt], 
                            BIC=[BIC_pkfit_2cmp_comb, BIC_pkfit_base_wt], 
                            OFV=[minustwoLL_pkfit_2cmp_comb, minustwoLL_pkfit_base_wt])
#


## post-processing 


## Save coefficients for model parameter estimates to a text document 
open("WCOP//poppk_data//base_wt//coef_base_wt.txt", "w") do io
  println(io, coeftable(pkfit_base_wt))
end


## Make an inference
infer_base_wt = infer(pkfit_base_wt) |> coeftable
CSV.write("WCOP//poppk_data//base_wt//infer_base_wt.csv", infer_base_wt)

## Generate diagnostics
inspect_base_wt = inspect(pkfit_base_wt) |> DataFrame
CSV.write("WCOP//poppk_data//base_wt//inspect_base_wt.csv", inspect_base_wt)

## post-hoc concentration predictions 
pred_base_wt = DataFrame(predict(pkfit_base_wt, obstimes=0:0.5:168))
CSV.write("WCOP//poppk_data//base_wt//pred_base_wt.csv", pred_base_wt)

## individual predictions
indiv_base_wt = reduce(vcat, DataFrame.(icoef(pkfit_base_wt)))
CSV.write("WCOP//poppk_data//base_wt//indiv_base_wt.csv", indiv_base_wt)

## empirical bayes 
ebes_base_wt = DataFrame(empirical_bayes(pkfit_base_wt))
CSV.write("WCOP//poppk_data//base_wt//ebes_base_wt.csv", ebes_base_wt)

## residuals 
resid_base_wt = DataFrame(wresiduals(pkfit_base_wt))
CSV.write("WCOP//poppk_data//base_wt//resid_base_wt.csv", resid_base_wt)

## shrinkage 
open("WCOP//poppk_data//base_wt//shrinkage_base_wt.txt", "w") do io
  println(io, ηshrinkage(pkfit_base_wt))
end
















###################### CRCL ##########################################


mdl_base_wt_crcl = @model begin
    @param begin
      tvcl   ∈ RealDomain(lower = 0.0001)
      tvvc   ∈ RealDomain(lower = 0.0001)
      tvq    ∈ RealDomain(lower = 0.0001)
      tvvp   ∈ RealDomain(lower = 0.0001)
      Ω      ∈ PDiagDomain(3)
      σ_add  ∈ RealDomain(lower = 0.0001)
      σ_prop ∈ RealDomain(lower = 0.0001)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @covariates WT CRCL 

    @pre begin
        wtCL = (WT/70)^0.75
        wtV  = (WT/70)
        crcl_eff = (CRCL/95)^0.75
        CL   = tvcl * crcl_eff * wtCL * exp(η[1])
        Vc   = tvvc * wtV  * exp(η[2])
        Q    = tvq  * wtCL * exp(η[3])
        Vp   = tvvp * wtV  
    end

    @dynamics Central1Periph1

    @derived begin
      CONC := @. Central/Vc
      DV ~ @. Normal(CONC, (abs(CONC)*σ_prop)+σ_add)
    end
end



params_base_wt_crcl = (
    tvvc   = 5,
    tvcl   = 0.05,
    tvq    = 0.05,
    tvvp   = 10,
    Ω      = Diagonal([0.01, 0.01, 0.01]),
    σ_add  = 0.1,
    σ_prop = 0.1)
#

## Maximum likelihood estimation
pkfit_base_wt_crcl = fit(mdl_base_wt_crcl,
              pop,
              params_base_wt_crcl,
              Pumas.FOCEI(),
              ensemblealg = EnsembleThreads())
#

serialize("WCOP//poppk_data//base_wt_crcl//pkfit_base_wt_crcl.jls", pkfit_base_wt_crcl)
#pkfit_base_wt_crcl = deserialize("./templates/data/base_wt_crcl/pkfit_base_wt_crcl.jls")


## model diagnstic 
LL_pkfit_base_wt_crcl = Pumas.marginal_nll(pkfit_base_wt_crcl)
minustwoLL_pkfit_base_wt_crcl = 2 * LL_pkfit_base_wt_crcl
AIC_pkfit_base_wt_crcl = aic(pkfit_base_wt_crcl)
BIC_pkfit_base_wt_crcl = bic(pkfit_base_wt_crcl) 


df_comp_base_wt_crcl = DataFrame(Model=["Base", "Base+WT", "Base+WT+CRCL"], 
                            AIC=[AIC_pkfit_2cmp_comb, AIC_pkfit_base_wt, AIC_pkfit_base_wt_crcl], 
                            BIC=[BIC_pkfit_2cmp_comb, BIC_pkfit_base_wt, BIC_pkfit_base_wt_crcl], 
                            OFV=[minustwoLL_pkfit_2cmp_comb, minustwoLL_pkfit_base_wt, minustwoLL_pkfit_base_wt_crcl])
#


## post-processing 


## Save coefficients for model parameter estimates to a text document 
open("WCOP//poppk_data//base_wt_crcl//coef_base_wt_crcl.txt", "w") do io
  println(io, coeftable(pkfit_base_wt_crcl))
end


## Make an inference
infer_base_wt_crcl = infer(pkfit_base_wt_crcl) |> coeftable
CSV.write("WCOP//poppk_data//base_wt_crcl//infer_base_wt_crcl.csv", infer_base_wt_crcl)

## Generate diagnostics
inspect_base_wt_crcl = inspect(pkfit_base_wt_crcl) |> DataFrame
CSV.write("WCOP//poppk_data//base_wt_crcl//inspect_base_wt_crcl.csv", inspect_base_wt_crcl)

## post-hoc concentration predictions 
pred_base_wt_crcl = DataFrame(predict(pkfit_base_wt_crcl, obstimes=0:0.5:168))
CSV.write("WCOP//poppk_data//base_wt_crcl//pred_base_wt_crcl.csv", pred_base_wt_crcl)

## individual predictions
indiv_base_wt_crcl = reduce(vcat, DataFrame.(icoef(pkfit_base_wt_crcl)))
CSV.write("WCOP//poppk_data//base_wt_crcl//indiv_base_wt_crcl.csv", indiv_base_wt_crcl)

## empirical bayes 
ebes_base_wt_crcl = DataFrame(empirical_bayes(pkfit_base_wt_crcl))
CSV.write("WCOP//poppk_data//base_wt_crcl//ebes_base_wt_crcl.csv", ebes_base_wt_crcl)

## residuals 
resid_base_wt_crcl = DataFrame(wresiduals(pkfit_base_wt_crcl))
CSV.write("WCOP//poppk_data//base_wt_crcl//resid_base_wt_crcl.csv", resid_base_wt_crcl)

## shrinkage 
open("WCOP//poppk_data//base_wt_crcl//shrinkage_base_wt_crcl.txt", "w") do io
  println(io, ηshrinkage(pkfit_base_wt_crcl))
end

















##########################################################
#                                                        #
#             FINAL MODEL QUALIFICATION                  #
#                                                        #
##########################################################