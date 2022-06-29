using Pumas
using CSV
using DataFramesMeta
using PumasUtilities
using CairoMakie
using AlgebraOfGraphics

# 1. Read the dataset in
pkpd_df = CSV.read("pkpd_data/pkpd_idr1.csv", DataFrame)

# 2. Map to read_pumas for population. Only use PK
pkpdpop = read_pumas(pkpd_df, observations = [:dv], covariates = [:BSL])

# 3. Plot Data
@chain pkpd_df begin
  dropmissing(:dv)
  data(_) *
  mapping(:time, :dv, layout = :id => nonnumeric) *
  visual(ScatterLines, linewidth = 4, markersize = 8, markercolor = :blue)
  draw(axis = (xlabel = "Time (hours)", ylabel = "PK Concentration", yscale = log10),
       figure = (resolution = (1400, 1000), fontsize = 22))
end

# 4. Specify pk model and parameters
inf_2cmt = @model begin
    @param   begin
        tvcl ∈ RealDomain(lower=0)
        tvvc ∈ RealDomain(lower=0)
        tvq ∈ RealDomain(lower=0)
        tvvp ∈ RealDomain(lower=0)
        Ω_pk ∈ PDiagDomain(4)
        σ_prop_pk ∈ RealDomain(lower=0)
    end

    @random begin
      ηpk ~ MvNormal(Ω_pk)
    end

    @pre begin
        CL = tvcl * exp(ηpk[1])
        Vc = tvvc * exp(ηpk[2])
        Q = tvq   * exp(ηpk[3])
        Vp = tvvp * exp(ηpk[4])
    end

    @dynamics begin
       Central' =  - (CL/Vc)*Central + (Q/Vp)*Peripheral - (Q/Vc)*Central
       Peripheral' = (Q/Vc)*Central - (Q/Vp)*Peripheral
     end

    @derived begin
        cp := @. (Central/Vc)
        dv ~ @. Normal(cp, sqrt(cp^2*σ_prop_pk))
    end
end

## specify the parameters
twocomp_params = (tvcl = 1.5,
                tvvc = 25.0,
                tvq = 5.0,
                tvvp = 150.0,
                Ω_pk = Diagonal([0.05,0.05,0.05,0.05]),
                σ_prop_pk = 0.02)
#

# 5. Fit the data
pkfit = fit(inf_2cmt,
            pkpdpop,
            twocomp_params,
            Pumas.FOCE())

# 6. explore the results
pkinspect = inspect(pkfit)
goodness_of_fit(pkinspect, ols = false)

# Assume that the pk model above is the final model

# 7. Lets extract the individual parameters
indpars = DataFrame(pkinspect.icoefs)
select!(indpars, :id, :CL, :Vc, :Q, :Vp)
indpars[!, :id] = parse.(Int64, indpars.id)

# 8. Merge the individual parameters with original dataset
pkparam_pkpddf = @chain pkpd_df begin
  leftjoin(indpars, on = [:id])
  dropmissing([:CL, :Vc, :Q, :Vp])
  rename(:CL => :iCL, :Vc => :iVc, :Q => :iQ, :Vp => :iVp)
end

# 9. Convert data into population
pdpop = read_pumas(pkparam_pkpddf, observations = [:resp], 
                  covariates = [:BSL, :iCL, :iVc, :iQ, :iVp])

# 10. Specify the PKPD model
inf_2cmt_lin_turnover = @model begin
    @param   begin
        # PD parameters
        tvturn ∈ RealDomain(lower=0)
        tvebase ∈ RealDomain(lower=0)
        tvec50 ∈ RealDomain(lower=0)
        Ω_pd ∈ PDiagDomain(1)
        σ_add_pd ∈ RealDomain(lower=0)

    end

    @random begin
      ηpd ~ MvNormal(Ω_pd)
    end

    @covariates iCL iVc iQ iVp

    @pre begin
        CL = iCL
        Vc = iVc
        Q  = iQ
        Vp = iVp

        ebase = tvebase*exp(ηpd[1])
        ec50 = tvec50
        emax = 1
        turn = tvturn
        kout = 1/turn
        kin0 = ebase*kout
    end

    @init begin
        Resp = ebase
    end

    @vars begin
        conc := Central/Vc
        edrug := emax*conc/(ec50 + conc)
        kin := kin0*(1-edrug)
    end

     @dynamics begin
         Central' =  - (CL/Vc)*Central + (Q/Vp)*Peripheral - (Q/Vc)*Central
         Peripheral' = (Q/Vc)*Central - (Q/Vp)*Peripheral
         Resp' = kin - kout*Resp
     end

    @derived begin
        resp ~ @. Normal(Resp, sqrt(σ_add_pd))
    end
end
## specify the parameters
turnover_params = (tvturn = 10,
                tvebase = 10,
                tvec50 = 0.3,
                Ω_pd = Diagonal([0.05]),
                σ_add_pd = 0.2)

# 11. Fit the PD data via sequential pkpd
pkpdfit = fit(inf_2cmt_lin_turnover,
            pdpop,
            turnover_params,
            Pumas.FOCE())

# 12. explore the results
pdinspect = inspect(pkpdfit)
goodness_of_fit(pdinspect, ols = false)

