using CSV, Pumas, PumasUtilities, Bioequivalence, Dates

# Read Data
pkdata = CSV.read("nca_data/data4NCA_sad.csv", DataFrame; missingstring="")

#  map NCA data from CSV
ncapop = read_nca(pkdata; observations=:concentration, group=[:doselevel])

# plot means
sp = summary_observations_vs_time(
    ncapop;
    paginate=true,
    columns=3,
    rows=1,
    axis=(; xlabel="Time (hours)", ylabel="Concentration (μg/mL)"),
    facet=(; combinelabels=true),
    figure=(; resolution=(1400, 800), fontsize=25),
)
sp[1]

ot = observations_vs_time(
    ncapop;
    axis=(;
        xlabel="Time (hours)",
        ylabel="CTMX Concentration (μg/mL)",
        yscale=log10,
        ytickformat=x -> string.(round.(x; digits=1)),
        yticks=[0.1, 10, 20, 30, 50, 100, 1000],
    ),
    columns=4,
    rows=3,
    facet=(; combinelabels=true),
)

# NCA for specific parameters

# Clearance
NCA.cl(ncapop)

# Cmax
NCA.cmax(ncapop)

# AUC0-2
NCA.auc(ncapop; interval=(0, 2))

# lambdaz
NCA.lambdaz(ncapop; interval=(48, 72))

# Generating report for NCA

# run the complete NCA analysis
pk_nca = run_nca(
    ncapop;
    sigdigits=3,
    studyid="STUDY-001",
    studytitle="SAD",
    author=[("Auth1", "Company1"), ("Auth2", "Company2")],
    sponsor="Company",
    date=Dates.now(),
    conclabel="Concentration (μg/mL)",
    grouplabel="Dose (mg)",
    timelabel="Time (Hr)",
    versionnumber=v"0.1",
    header="Pumas-AI",
    footer="Confidential",
)

# Plots

# Lambdaz fits
sf = subject_fits(
    ncapop;
    paginate=true,
    axis=(; xlabel="Time (hours)", ylabel="Concentration (μg/mL)", yscale=log10),
    facet=(; combinelabels=true),
)
sf[1]

# parameter distribution plots

# cl_f_obs distribution plots
parameters_vs_group(pk_nca; parameter=:cl_f_obs)

# cmax_f_obs distribution plots
parameters_vs_group(pk_nca; parameter=:cmax)
