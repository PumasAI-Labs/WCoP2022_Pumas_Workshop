using Random
using Pumas
using PumasUtilities
using CairoMakie
using DataFramesMeta

###########################
# Customizing Pumas Plots #
###########################

########################
## First order aborption

# Create subjects - dosing regimen + covariates
dose = DosageRegimen(100; time=0) # By default, cmt=1 (https://docs.pumas.ai/stable/basics/doses_subjects_populations/#Dosage-Regimen-Terminology)
choose_covariates() = (;
                       wt=rand(Random.seed!(1234), 55:80),
                       dose = 100
)
subj_with_covariates = map(1:10) do i
    Subject(
            id=i,
            events=dose,
            covariates=choose_covariates(),
            observations=(; conc = nothing)
    )
end

# Create model
FOabs = @model begin
    @param begin
        tvcl ∈ RealDomain(; lower=0)
        tvvc ∈ RealDomain(; lower=0)
        tvka ∈ RealDomain(; lower=0)
        Ω    ∈ PDiagDomain(3)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates wt

    @pre begin
        CL = tvcl * (wt / 70)^0.75 * exp(η[1])
        Vc = tvvc * (wt / 70) * exp(η[2])
        Ka = tvka * exp(η[3])
    end

    @dynamics Depots1Central1

    @observed begin
        conc = @. Central/Vc
    end
end

# Model parameters for simulation
param = (;
         tvcl=5, tvvc=20, tvka=1,
         Ω=Diagonal([0.04, 0.04, 0.04])
)

# Perform simulations         
sims = simobs(FOabs, subj_with_covariates, param; obstimes=0:.1:24, rng=Random.seed!(123))

# Simulations can be converted to a tabular format for exporting
sims_df = DataFrame(sims)

# Plot simulated profiles using sim_plot
# Let's customize this?
# Axis: https://makie.juliaplots.org/stable/documentation/api_reference/#Axis
sim_plot(sims; observations=[:conc], 
         figure=(; fontsize = 18), 
         axis=(;
               xlabel = "Time (hr)", 
               ylabel = "Predicted Concentration (ng/mL)", 
               title = "First-order absorption")
)    

# More Pumas `sim_plot` custmization:
# https://docs.pumas.ai/stable/analysis/plots/#PumasPlots.sim_plot

###############################
#   Layouting and Subplots    #
###############################

# We can do layouts with a Figure() from Makie
f = Figure() # empty Figure
# Multiple dispatch on Makie's plotting functions
# 🤯
scatter(f[1, 1], 1:10)           # plotting on row 1 col 1
lines(f[1, 2], 1:10)             # plotting on row 1 col 2
heatmap(f[2, 1:2], rand(10, 10)) # plotting on row 2 cols 1 to 2

# We can do that with *ANY* Pumas plot:
# https://docs.pumas.ai/stable/analysis/plots/#Customizing-Plots

