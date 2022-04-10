#verbs (read_pumas, fit, infer, inspect, evaluate_diagnostics)
#nouns (Subject, Population, DosageRegimen)


m1 = @model begin
    @metadata begin
        desc = "base model: 1comp"
        timeu = u"hr"
    end

    @param begin
        "Clearance (L/hr)"
        tvcl ∈ RealDomain(lower=0.0001)
        "Volume (L)"
        tvvc ∈ RealDomain(lower=0.0001)
        """
        - ΩCL
        - ΩVc
        """
        Ω ∈ PDiagDomain(2) #covariance matrix 
        # σ_add  ∈ RealDomain(lower = 0.0001) # standard deviation 
        # σ_prop ∈ RealDomain(lower = 0.0001) # standard deviation 
        "Additive RUV"
        σ²_add ∈ RealDomain(lower=0.0001) # variance 
        "Proportional RUV"
        σ²_prop ∈ RealDomain(lower=0.0001) # variance  
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates crcl sex

    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvvc * exp(η[2])
    end

    @dynamics Central1

    @derived begin
        CONC = @. Central / Vc #suppressing CONC from being output using ":="
        """
        DrugY Concentration (ng/mL)
        """
        DV ~ @. Normal(CONC, √(abs(CONC)^2 * σ²_prop + σ²_add)) # using variance
    end
end



params_1cmp_comb = (
    tvvc=5,
    tvcl=0.2,
    Ω=Diagonal([0.01, 0.01]),
    σ²_add=0.01,
    σ²_prop=0.01)
#

s1 = Subject(id="Allison",
    covariates=(crcl=90, sex="f"),
    events=DosageRegimen(100, time=0, cmt=:Central, addl=5, ii=24))
s2 = Subject(id="Rahul",
    covariates=(crcl=150, sex="m"),
    events=DosageRegimen([200, 100], time=[0, 2], cmt=[:Central, :Central], addl=[0, 5], ii=[0, 24]))

pop = [s1, s2]

# simulate observations (simobs)
sims = simobs(m1, pop, params_1cmp_comb, obstimes=0:1:168)
sim_plot(sims, observations=[:CONC])

# population

a = [4, 5, 6]

log10.(a)

map(x -> log10(x), a)

pop2 = map(x -> Subject(id=x, events=DosageRegimen(100, addl=5, ii=24), covariates=(crcl=90, sex="f")), 1:10)
pop3 = map(x -> Subject(id=x, events=DosageRegimen(200, addl=3, ii=48), covariates=(crcl=90, sex="f")), 100:110)
pop = vcat(pop2, pop3)
# simulate observations (simobs)
sims = simobs(m1, pop, params_1cmp_comb, obstimes=0:1:168)
sim_plot(sims, observations=[:CONC])

# 100 mg PO qd X 3 days
ev1 = DosageRegimen(100, cmt=:Depot, time=0, addl=2, ii=24)
ev2 = DosageRegimen(50, cmt=:Depot,  time=0, addl=3, ii=12)
ev3 = DosageRegimen(ev1, ev2, offset = 3*24)

ev4 = DosageRegimen([100, 50], cmt=1, time=[0, 72], addl=[2, 3], ii=[24, 12])
s4 = Subject(id="Jose", events=ev4, covariates=(crcl=80, sex="m"))

sims = simobs(m1, s4, params_1cmp_comb, obstimes=0:1:168)
sim_plot(sims, observations=[:CONC])

simdf = DataFrame(sims)
simpop = read_pumas(simdf,
                    observations = [:DV])

# 
simpop2 = Subject(sims)                    