using Pumas
using CSV
using Random
using PumasUtilities
using CairoMakie


# Overview of model
# structural model - parameters (structure of cmt - defined by equations)
    # Ex: 1 cmt model with 1st order absorption -- have param of Ka, CL & V
    # can be differential equations or analytical solutions
# covariate model - relation between covariates & parameters
    # Ex: WT on CL (thus wt is a covarite of the model)
    # can be constant, time varying or it can be on an occasion
# variability model - between subject or within variability
    # variability of parameters in the population - have a distribution assigned to them
    # η ~ Normal(Θ, Ω) -- eta relate parameter to pop values
# error model - additive, proportional, combination, t-dist, etc.



############################## MODEL EXAMPLE ##################################


model = @model begin
# note - parameters you pass into model MUST be type NamedTuple
    @param begin #equal to the theta, sigma & omega block all together in Phoenix
        # EXAMPLE OF PROVIDING INTIAL ESTIMATE IN PARAM BLOCK: 
        #tvKa ∈ RealDomain(lower=0.0) # ∈ stands for "is in", init = 1 means initial estimate = 1
        tvCL ∈ RealDomain(lower=0.0) # real domain assigns a space for parameters to search during optimization
        tvVc ∈ RealDomain(lower=0.0) # can also pass initial estimates outside the model
        # these are the population parameters ^ (FIXED EFFECTS) - need to relate fixed effects to ranodm effects (BSV)
        # they must be linked to the BSV distribution (etas)
        # EXAMPLE OF PROVIDING INTIAL ESTIMATE IN PARAM BLOCK: 
        #tvKa ∈ RealDomain(lower=0.0, init=1.0) 
        #tvCL ∈ RealDomain(lower=0.0, init=4.0) 
        #tvVc ∈ RealDomain(lower=0.0, init=70.0) 
        Ω ∈ PDiagDomain(2) # #omega is a parameter to be estimated for BSV calculation of eta (this is specifying a diagnal matrix)
        #Ω ∈ PDiagDomain(init=[0.1,0.1,0.1]) # example with init 
        # Ω ∈ PSDDomain(2) specifies a correlation materix
        # σ_prop ∈ RealDomain(lower=0)
        # specify parameters in the model & there boundaries
        # EXAMPLE: want a correlation with CL and Vc, but not with Ka
            # Ωka ∈ PDiagDomain(1)
            # Ωclv ∈ PSDDomain(2)
        σ_prop ∈ RealDomain(lower=0.0)
    end

    @random begin
        η ~ MvNormal(Ω) # multivariate normal distribution with variance of omega
        # if you know your eta are heavily tailed to the end, can do other distributions
        # saying that eta came from distribution of omega which has 2 elements
        # thus eta will also have 2 elements (eta1 & eta2)
        # DIFFERENCE BETWEEN MVNORMAL & NORMA DISTRIBUTION:
        # EXAMPLE: want a correlation with CL and Vc, but not with Ka
            # ηka  ~ Normal(Θ, Ωka)
            # ηclv ~ MvNormal(Ωclv)
        # i.e. multivariate normal is when there is a correlation
    end

    @covariates WT #specifying the name of the covariate

    @pre begin
        #Ka = tvKa*exp(η[1]) # if oral -- would need to re-number etas 
        CL = tvCL*(WT/70)^0.75*exp(η[1]) #eta 1 is eta CL which is the 1st element of diagnal matrix
        Vc = tvVc*(WT/70)*exp(η[2])
        # this is what we write in the PK block of Phoenix
    end

    #@dynamics Depots1Central1 is equivalent to the statement below - 1st order absoprtion 1 cmt model
    # ^ if using this MUST use case sensitive as specified when you run Depots1Central1() in julia terminal
    @dynamics begin
        Central' =  -(CL/Vc)*Central #central compartment
        # oral route example 
        #Depot'   = -Ka*Depot # absorption compartment (left side is rate, right side should all be rates)
        #Central' =  Ka*Depot - (CL/Vc)*Central #central compartment
        # have 3 parameters - Ka, CL & VC - define these in the param block
    end

    @derived begin
        cp = @. (Central/Vc)
        dv ~ @. Normal(cp, abs(cp)*σ_prop) # mean of Cp and a variance of Cp*sigma (proportional error model)
        # Y = F+ F*esp  in nonmem is same as above
        # dv ~ @. Normal(cp, σ_add) # additive error
        # dv ~ @. Normal(cp, (cp*σ_prop + σ_add)) # combination error
        # dv ~ @. Bernoulli(logistic(logisticparam)) # can call all different standard models here - just specify the mean & variance
        # note: ~ is saying it comes from a distribution (it belongs to)
            # using = is for assignment, ~ is for distribution
    end
end










###############################################################################
#                       BUILDING A DOSAGE REGIMEN                             #
###############################################################################

# create a dosing event
# dose amount = 100
# 1st dose given at time = 0
# give 3 additional doses (thus 4 in total)
# give those additional in 24 hour intervals
ev = DosageRegimen(1000, time=0, addl=3, ii=4)

# create a subject that we are giving this dosage regimen to
# id = identifier can be numeric or character
# events --> pass in what events the subject will have (ev) - i.e. 4 doses in this scenario
# covariates --> can specify whatever covariates you will read into model 
s1 = Subject(id=1, events=ev, covariates=(WT=70,), observations=(dv=nothing,)) # note the "," after WT=70 (example of tricking the system)

# set parameters
#param = init_params(model) #would use this if chose to use "init" in pre block 
param = (tvCL = 1,
        tvVc = 10,
        Ω = Diagonal([0.1, 0.1]),
        σ_prop = 0.1)
#



# simulate every 1 hour until time 120 hours
# simulating for subject 1 with parameters above and the model previously described
subjobs = simobs(model, s1, param, obstimes=0:1:120)
sim_df = DataFrame(subjobs) # can convert to a dataframe if wish to save as CSV


# simple plot 
sim_plot(subjobs)



# customize (similar to ggplots)
sim_plot(model,
            [subjobs], observations = [:cp], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (ng/mL)"))
#








##############################################################################
#               CREATE MANY SUBJECTS WITH RANDOM COVARIATES                  #
##############################################################################

# function created to give random covariates in a simulated population
choose_covariates() = (WT = rand(55:80),)

# type of population - also a population when you combine subjects such as [s1, s2]
pop_with_covariates = map(i->Subject(id=i, events=ev,observations=(dv=nothing,),covariates=choose_covariates()),1:50)
# assuming all patients are getting the same event (ev)
# given random covariates from choose_covariate function



obs_pk = simobs(model, pop_with_covariates, param, obstimes=0:1:72)
# simulation output is readily available to refit for an observation routine
# i.e. have all the required elements for Pumas model
simdf = DataFrame(obs_pk) # type dataframe (not type Population)


# in pumas for any analysis, our dataset must be of type Subject or type population
# thus need to map this dataframe to a population type
# now that we have the dataframe of our simulated data, we can run our model


sim_plot(model,
            obs_pk, observations = [:cp], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (ng/mL)"))
#


sim_plot(model,
            obs_pk, observations = [:dv], 
            figure = (fontsize = 18, ), 
            axis = (xlabel = "Time (hr)", 
                    ylabel = "Predicted Concentration (ng/mL)"))
#

# example - pkdata = read_pumas(simdf, dvs=[:parent, :metabolite], cvs=[:isPM, :wt])


# fit -- fitting the data (same syntax as simobs)
pkres = fit(model, 
            Subject.(obs_pk), 
            param, 
            Pumas.FOCE(),
            ensemblealg = EnsembleThreads()) # type is FittedPumasModel
#
# adding "ensemblealg=EnsembleThreads()" is not a necessary component 
# 
# to view model estimates in REPL 
coef(pkres)

# to generate a dataframe of model estimates
coeftable(pkres)

# to make an inference on our fit - use infer (produces confidence intervals & SE)
# evaluating if the parameters we have are correct or not in terms of precision
pkinfer = infer(pkres)


# to get individual predicitions - use inspect 
insp = inspect(pkres)
insp_df = DataFrame(insp)

# to save any CSV file: 
CSV.write("inspect_pk.csv", insp_df)

