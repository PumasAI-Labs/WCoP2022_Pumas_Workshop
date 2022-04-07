# Developing a population PK model for iv dosing with covariates using Pumas

# Load the necessary Libraries
using CSV, Pumas, PumasUtilities, Bioequivalence, Dates

# Read data
pkdata = CSV.read("poppk_data/data4poppk_iv.csv", DataFrame; missingstring=".")

# Coverting the DataFrame to a collection of Subjects (Population)
population = read_pumas(
    pkdata;
    id=:ID,
    time=:TIME,
    observations=[:CONC],
    evid=:EVID,
    amt=:AMT,
    cmt=:CMT,
    covariates=[:AGE, :WEIGHT, :ISMALE, :eGFR],
)

# Model definition
model = @model begin

    @param begin
        # here we define the parameters of the model
        tvcl ∈ RealDomain(; lower=0.1) # typical clearance
        tvvc ∈ RealDomain(; lower=1.0) # typical central volume of distribution
        Ω ∈ PDiagDomain(2)             # between-subhect variability
        σ ∈ RealDomain(; lower=0.0)    # residual error
    end

    @random begin
        # here we define random effects
        η ~ MvNormal(Ω) # multi-variate Normal with mean 0 and covariance matrix Ω
    end

    @pre begin
        # pre computations and other statistical transformations
        CL = tvcl * exp(η[1])
        Vc = tvvc * exp(η[2])
    end

    # here we define compartmends and dynamics
    @dynamics Central1 # same as Central' = -(CL/Vc)*Central (see Pumas documentation)

    @derived begin
        # here is where we calculate concentration and add residual variability
        # tilde (~) means "distributed as"
        cp = @. Central / Vc
        CONC ~ @. Normal(cp, sqrt(σ))
    end
end

# Parameter values
params = (tvcl=1.0, tvvc=10.0, Ω=Diagonal([0.09, 0.09]), σ=3.16)

# Fit covariate model
fit_results = fit(model, population, params, Pumas.FOCE())
fit_infer = infer(fit_results)

fit_inspect = inspect(fit_results)
fit_diagnostics = evaluate_diagnostics(fit_inspect)

fit_vpc = vpc(fit_results)

fit_diagnostics = evaluate_diagnostics((fit_inspect, fit_vpc),)
