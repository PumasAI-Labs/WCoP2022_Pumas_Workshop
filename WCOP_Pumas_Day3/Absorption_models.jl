
##################
## Load Libraries

using Random
using Pumas
using PumasUtilities
using CairoMakie
using DataFramesMeta

Random.seed!(1234)  # Set random seed for reproducibility

########################
## First order aborption

# Create subjects - dosing regimen + covariates
dose = DosageRegimen(100, time = 0) # By default, cmt=1 (https://docs.pumas.ai/stable/basics/doses_subjects_populations/#Dosage-Regimen-Terminology)
choose_covariates() = (wt = rand(55:80), dose = 100)
subj_with_covariates = map(1:10) do i
    Subject(id = i,
            events = dose,
            covariates = choose_covariates(),
            observations = (conc = nothing,))
end

# Create model
FOabs = @model begin
    @param begin
        tvcl ∈ RealDomain(lower=0)
        tvvc ∈ RealDomain(lower=0)
        tvka ∈ RealDomain(lower=0)
        Ω    ∈ PDiagDomain(3)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates wt

    @pre begin
        CL = tvcl*(wt/70)^0.75*exp(η[1])
        Vc = tvvc*(wt/70)*exp(η[2])
        Ka = tvka*exp(η[3])
    end

    @dynamics Depots1Central1

    @observed begin
        conc = @. Central/Vc
    end
end

# Model parameters for simulation
param = (tvcl = 5, tvvc = 20, tvka = 1,
         Ω = Diagonal([0.04, 0.04, 0.04]))

# Perform simulations         
sims = simobs(FOabs, subj_with_covariates, param, obstimes = 0:.1:24)

# Simulations can be converted to a tabular format for exporting
sims_df = DataFrame(sims)

# Plot simulated profiles using sim_plot
sim_plot(sims, observations =[:conc], 
        figure = (fontsize = 18, ), 
        axis = (xlabel = "Time (hr)", 
                ylabel = "Predicted Concentration (ng/mL)", 
                title = "First-order absorption"))         

########################
## Zero order aborption                

# Create subjects - dosing regimen + covariates
dose = DosageRegimen(100, time = 0, rate = -2)  # Note: rate must be -2 for duration modeling

subj_with_covariates = map(1:10) do i
    Subject(id = i,
            events = dose,
            covariates = choose_covariates(),
            observations = (conc = nothing,))
end

# Create model
ZOabs = @model begin
    @param begin
        tvcl  ∈ RealDomain(lower=0)
        tvvc  ∈ RealDomain(lower=0)
        tvdur ∈ RealDomain(lower=0)
        Ω     ∈ PDiagDomain(3)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates wt

    @pre begin
        CL = tvcl*(wt/70)^0.75*exp(η[1])
        Vc = tvvc*(wt/70)*exp(η[2])
    end

    @dosecontrol begin
        duration = (Central = tvdur*exp(η[3]),)
    end

    @dynamics Central1

    @observed begin
        conc = @. Central/Vc
    end
end

# Model parameters for simulation
param = (tvcl = 0.792, tvvc = 13.7, tvdur = 5.0,
         Ω = Diagonal([0.04, 0.04, 0.04]))

# Perform simulations         
sims = simobs(ZOabs, subj_with_covariates, param, obstimes = 0:0.1:48)

# Convert to dataframe 
sims_df = DataFrame(sims)

# Plot simulated profiles
sim_plot(sims, observations =[:conc], 
        figure = (fontsize = 18, ), 
        axis = (xlabel = "Time (hr)", 
                ylabel = "Predicted Concentration (ng/mL)", 
                title = "Zero-order absorption"))         

###############################################
## Parallel First Order + Zero order absorption 

# Create model
ZOFO_paral_abs = @model begin

    @param begin
        tvcl  ∈ RealDomain(lower=0)
        tvvc  ∈ RealDomain(lower=0)
        tvka  ∈ RealDomain(lower=0)
        tvdur ∈ RealDomain(lower=0)
        tvbio ∈ RealDomain(lower=0)
        tvlag ∈ RealDomain(lower=0)
        Ω     ∈ PDiagDomain(2)
    end
  
    @random begin
        η ~ MvNormal(Ω)
    end
  
    @covariates wt
  
    @pre begin
        CL = tvcl*(wt/70)^0.75*exp(η[1])
        Vc = tvvc*(wt/70)*exp(η[2])
        Ka = tvka
    end
  
    @dosecontrol begin
        duration = (Central = tvdur,)
        bioav    = (Depot = tvbio, Central = 1 - tvbio)
        lags     = (Depot = tvlag,)
    end
  
    @dynamics begin
        Depot'   = -Ka*Depot
        Central' =  Ka*Depot - (CL/Vc)*Central
    end
  
    @observed begin
        conc = @. Central/Vc
    end
  
end

# Model parameters for simulation
param = (tvcl = 5, tvvc = 50, tvka = 1.2, 
         tvdur = 2, tvbio = 0.5, tvlag = 1,
         Ω = Diagonal([0.04, 0.04]))

# Create subjects - dosing regimen + covariates
dose_zo = DosageRegimen(100, time = 0, cmt = 2, rate = -2, evid = 1)  # Note: rate must be -2 for duration modeling
dose_fo = DosageRegimen(100, time = 0, cmt = 1, rate = 0,  evid = 1)
dose    = DosageRegimen(dose_fo, dose_zo)  # Actual dose is made up of 2 virtual doses

dose    = DosageRegimen([100, 100], time = 0, cmt = [1, 2], rate=[0, -2], evid=1) # Another way to create the above dosing regimen

subj_with_covariates = map(1:10) do i
    Subject(id = i,
            events = dose,
            covariates = choose_covariates(),
            observations = (conc = nothing,))
end

# Perform simulations         
sims = simobs(ZOFO_paral_abs, subj_with_covariates, param, obstimes = 0:0.1:24)

# Plot simulated profiles
sim_plot(sims, observations =[:conc], 
        figure = (fontsize = 18, ), 
        axis = (xlabel = "Time (hr)", 
                ylabel = "Predicted Concentration (ng/mL)", 
                title = "Zero- and first-order parallel absorption"))         

#####################################
## Two Parallel First Order Processes

# Create model
two_parallel_foabs = @model begin
    @param begin
        tvcl  ∈ RealDomain(lower=0)
        tvvc  ∈ RealDomain(lower=0)
        tvka1 ∈ RealDomain(lower=0)
        tvka2 ∈ RealDomain(lower=0)
        tvlag ∈ RealDomain(lower=0)
        tvbio ∈ RealDomain(lower=0)
        Ω     ∈ PDiagDomain(6)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates wt

    @pre begin
        CL  = tvcl*(wt/70)^0.75*exp(η[1])
        Vc  = tvvc*(wt/70)*exp(η[2])
        Ka1 = tvka1*exp(η[3])
        Ka2 = tvka2*exp(η[4])
    end

    # @dosecontrol begin
    #   lags  = (SR = tvlag*exp(η[5]),)
    #   bioav = (IR = tvbio*exp(η[6]), SR = (1 - tvbio)*exp(η[6]))
    # end

    @dosecontrol begin
        lags  = (Depot1 = tvlag*exp(η[5]),)
        bioav = (Depot2 = tvbio*exp(η[6]), Depot1 = (1 - tvbio)*exp(η[6]))
    end

    # @dynamics begin
    #   IR'      = -Ka1*IR
    #   SR'      = -Ka2*SR
    #   Central' =  Ka1*IR + Ka2*SR - Central*CL/Vc
    # end
    @dynamics Depots2Central1

    @observed begin
        conc = @. Central/Vc
    end
end

# Model parameters for simulation
param = (tvcl = 5, tvvc = 50, tvka1 = 0.8,
         tvka2 = 0.6, tvlag = 5, tvbio = 0.5,
         Ω = Diagonal([0.04, 0.04, 0.36, 0.36, 0.04, 0.04]))

# Create subjects - dosing regimen + covariates
dose_fo1 = DosageRegimen(100, cmt = 1, time = 0)
dose_fo2 = DosageRegimen(100, cmt = 2, time = 0)
dose     = DosageRegimen(dose_fo1, dose_fo2)  # Actual dose is made up of 2 virtual doses

subj_with_covariates = map(1:10) do i
    Subject(id = i,
            events = dose,
            covariates = choose_covariates(),
            observations = (conc = nothing,))
end

# Perform simulations  
sims = simobs(two_parallel_foabs, subj_with_covariates, param, obstimes = 0:.1:48)

# Plot simulated profiles
sim_plot(sims, observations =[:conc], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (ng/mL)", 
                    title = "Two Parallel first-order absorption"))

#####################
## Weibull absorption

# Create model
weibullabs = @model begin
    @param begin
        tvcl ∈ RealDomain(lower=0)
        tvvc ∈ RealDomain(lower=0)
        tvka ∈ RealDomain(lower=0)
        tvγ  ∈ RealDomain(lower=0)
        Ω    ∈ PDiagDomain(4)
    end
  
    @random begin
         η ~ MvNormal(Ω)
    end
  
    @covariates wt
  
    @pre begin
        CL  = tvcl*(wt/70)^0.75*exp(η[1])
        Vc  = tvvc*(wt/70)*exp(η[2])
        Ka∞ = tvka*exp(η[3])           # Maximum Ka as t → ∞
        γ   = tvγ*exp(η[4])            # Controls the steepness of the Ka curve
        Kaᵗ = 1 - exp(-(Ka∞*t)^γ)      # Weibull function
    end
  
    @dynamics begin
        Depot'   = -Kaᵗ*Depot
        Central' =  Kaᵗ*Depot - (CL/Vc)*Central
    end
  
    @derived begin
        conc = Central/Vc
    end
  
end

# Model parameters for simulation
param = (tvcl = 5, tvvc = 50, tvka = 0.4, tvγ = 4,
         Ω = Diagonal([0.04, 0.04, 0.36, 0.04]))

# Create subjects - dosing regimen + covariates
dose = DosageRegimen(100, cmt = 1, time = 0)

subj_with_covariates = map(1:10) do i
    Subject(id = i,
            events = dose,
            covariates = choose_covariates(),
            observations = (conc = nothing,))
end

# Perform simulations  
sims = simobs(weibullabs, subj_with_covariates, param, obstimes = 0:.1:24)

# Plot simulated profiles
sim_plot(sims, observations =[:conc], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (ng/mL)", 
                    title = "Weibull absorption"))

################################################
## Erlang absorption (Transit compartment model)

# Create model
erlangabs = @model begin
    @param begin
        tvcl  ∈ RealDomain(lower=0)
        tvvc  ∈ RealDomain(lower=0)
        tvktr ∈ RealDomain(lower=0)
        Ω     ∈ PSDDomain(3)
    end
  
    @random begin
          η ~ MvNormal(Ω)
    end
  
    @covariates wt
  
    @pre begin
            CL  = tvcl*(wt/70)^0.75*exp(η[1])
            Vc  = tvvc*(wt/70)*exp(η[2])
            Ktr = tvktr*exp(η[3])
    end
  
    @dynamics begin
        Depot'    = -Ktr*Depot
        Transit1' = Ktr*Depot    - Ktr*Transit1
        Transit2' = Ktr*Transit1 - Ktr*Transit2
        Transit3' = Ktr*Transit2 - Ktr*Transit3
        Transit4' = Ktr*Transit3 - Ktr*Transit4
        Transit5' = Ktr*Transit4 - Ktr*Transit5
        Central'  = Ktr*Transit5 - (CL/Vc)*Central
    end
  
    @observed begin
        conc = @. Central/Vc
    end
  
end

# Model parameters for simulation
param = (tvcl = 7, tvvc = 32, tvktr = 2.6,
         Ω = Diagonal([0.09, 0.22, 0.10, 0.52]))

# Create subjects - dosing regimen + covariates
dose = DosageRegimen(100, time = 0)

subj_with_covariates = map(1:10) do i
    Subject(id = i,
            events = dose,
            covariates = choose_covariates(),
            observations = (conc = nothing,))
end

# Perform Simulations
sims = simobs(erlangabs, subj_with_covariates, param, obstimes = 0:1:24)

# Plot simulated profiles
sim_plot(sims, observations =[:conc], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (ng/mL)", 
                    title = "Erlang absorption"))


##################################

#################
## Round tripping

# Create subjects - dosing regimen + covariates
dose = DosageRegimen(100, time = 0)
choose_covariates() = (wt = rand(55:80), dose = 100)
subj_with_covariates = map(1:10) do i
    Subject(id = i,
            events = dose,
            covariates = choose_covariates(),
            observations = (dv = nothing,))
end

# Create model
IVbolus = @model begin
    @param begin
        tvcl ∈ RealDomain(lower=0)
        tvvc ∈ RealDomain(lower=0)
        Ω    ∈ PDiagDomain(2) 
        σ_prop ∈ RealDomain(lower=0)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates wt

    @pre begin
        CL = tvcl*(wt/70)^0.75*exp(η[1])
        Vc = tvvc*(wt/70)*exp(η[2])
    end

    @dynamics Central1

    @derived begin
        conc := @. Central/Vc
        dv   ~ @. Normal(conc, abs(conc) * σ_prop)
    end
end

# Model parameters for simulation
param = (tvcl = 1, tvvc = 10, 
         Ω = Diagonal([0.09, 0.09]),
         σ_prop = 0.3)

# Perform simulations         
sims = simobs(IVbolus, subj_with_covariates, param, obstimes = 0:3:24)
sim_plot(sims, observations = [:dv])

# Convert simulations to a population 
pop = Subject.(sims)
# pop = read_pumas(DataFrame(sims), observations = [:dv], covariates=[:wt, :dose])

# Re-fit
param = (tvcl = 3, tvvc = 15, 
         Ω = Diagonal([0.1, 0.2]),
         σ_prop = 0.1)
refit_1cmp = fit(IVbolus, pop, param, Pumas.FOCE())

loglikelihood(IVbolus, pop, param, Pumas.FOCE())
findinfluential(IVbolus, pop, param, Pumas.FOCE())

######################
## Multiple trials

ntrials = 50

param = coef(refit_1cmp)
sims_ntrials   = map(rep -> simobs(IVbolus, subj_with_covariates, param, obstimes = 0:3:24), 1:ntrials)
#pops_ntrials   = map(sim -> read_pumas(DataFrame(sim), observations = [:dv], covariates = [:wt, :dose]), sims_ntrials)
pops_ntrials   = map(sim -> Subject.(sim), sims_ntrials)
refits_ntrials = map(pop -> fit(IVbolus, pop, param, Pumas.FOCE()), pops_ntrials)
infers_ntrials = map(fpm -> infer(fpm), refits_ntrials)

# Get all typical clearances
coeftable(refits_ntrials[1])
refits_clearances = map(refits_ntrials) do fpm
    df = coeftable(fpm)
    @rsubset! df :parameter == "tvcl"
    return df.estimate[1]
end

# Get CI for inference
lci_clearances = map(infers_ntrials) do fpmi
    df = coeftable(fpmi)
    @rsubset! df :parameter == "tvcl"
    df.ci_lower[1]/param.tvcl
end
power = sum(lci_clearances .> 0.6)/ntrials