using Pumas
# Create subjects - dosing regimen + covariates
dose = DosageRegimen(100, time = 0)
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
        σ_prop ∈ RealDomain(lower=0)
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

    @derived begin
        conc := @. Central/Vc
        dv   ~ @. Normal(conc, abs(conc) * σ_prop)
    end
end

# Model parameters for simulation
param = (tvcl = 5, tvvc = 20, tvka = 1,
         Ω = Diagonal([0.04, 0.04, 0.04]),
         σ_prop = 0.1)

# Perform simulations         
sims = simobs(FOabs, subj_with_covariates, param, obstimes = 0:.1:24)

# Convert simulations to a population 
pop = Subject.(sims, name=:dv)

# Re-fit
refit_1cmp = fit(FOabs, pop, param, Pumas.FOCE())