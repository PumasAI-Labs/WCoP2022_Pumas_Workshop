## SUMMARY OF NCA 
# 1. Read in the CSV file into a DataFrame 
# 2. Map DataFrame to type NCAPopulation -- read_nca() 
# 3. NCA Population use to run NCA -- run_nca() 
# 4. Summarize data of interest -- summarize()
# 5. Manually output results OR run a report -- report()


# Call Necessary Packages 
using CSV
using Chain
using DataFrames
using Dates
using NCA
using NCAUtilities
using NCA.Unitful
using PumasUtilities





###################################################################
#                                                                 #
#                  IV BOLUS SINGLE DOSE EXAMPLE                   #
#                                                                 #
###################################################################



# Load Data
df_bolus_sd = CSV.read("WCOP//nca_data//iv_bolus_sd.csv", DataFrame, missingstrings=["NA", ".", ""])



# Define Units 
timeu = u"hr"
concu = u"mg/L"
amtu  = u"mg"



# Map Dataframe to NCA Population
pop_bolus_sd = read_nca(df_bolus_sd,
                id            = :id,
                time          = :time,
                observations  = :conc,
                amt           = :amt,
                route         = :route,
                #timeu         = true, 
                #amtu          = true, 
                #concu         = true, 
                llq           = 0.001)
#




# Preview Data 
## individual plot - linear scale 
obsvstimes = observations_vs_time(pop_bolus_sd[1])

## individual plot - semi-log scale
obsvstimes = observations_vs_time(pop_bolus_sd[1], axis = (yscale = log,))

## grid of individual plots
ctplots = observations_vs_time(pop_bolus_sd[1:9], 
                                axis = (xlabel = "Time (hour)", 
                                        ylabel = "Drug Concentration (mg/L)"),
                                paginate = true, #creates multiple pages  
                                columns = 3, rows = 3, #number of col/rows per page 
                                facet = (combinelabels = true,)) #creates 1 label for each page
ctplots[1]

## mean concentration-time curve of population 
summary_observations_vs_time(pop_bolus_sd,
                                axis = (xlabel = "Time (hour)", 
                                ylabel = "Drug Concentration (mg/L)"))
#                        






# Perform Simple NCA
nca_bolus_sd = run_nca(pop_bolus_sd, sigdigits=3)





# Generate Individual NCA Parameters 
vz        = NCA.vz(pop_bolus_sd, sigdigits=3)  # Volume of Distribution/F, in this case since the drug is given orally
cl        = NCA.cl(pop_bolus_sd, sigdigits=3)  # Clearance/F, in this case since the drug is given orally
lambdaz   = NCA.lambdaz(pop_bolus_sd, threshold=3, sigdigits=3)  # Terminal Elimination Rate Constant, threshold=3 specifies the max no. of time point used for calculation
thalf     = NCA.thalf(pop_bolus_sd, sigdigits=3) # Half-life calculation
cmax_d    = NCA.cmax(pop_bolus_sd, normalize=true, sigdigits=3) # Dose Normalized Cmax
mrt       = NCA.mrt(pop_bolus_sd, sigdigits=3) # Mean residence time
aumc      = NCA.aumc(pop_bolus_sd, method=:linlog, sigdigits=3) # AUMC calculation, using :linlog method
individual_params      = innerjoin(vz,cl,lambdaz,thalf,cmax_d,mrt,aumc, on=[:id], makeunique=true)

# Other AUC calculation options 
auc0_12   = NCA.auc(pop_bolus_sd, interval=(0,12), method=:linuplogdown, sigdigits=3) #various other methods are :linear, :linlog
auc12_24  = NCA.auc(pop_bolus_sd, interval=(12,24), method=:linuplogdown, sigdigits=3) #looking at auc 12 to 24 hours (can make this interval anything!)
partial_aucs = NCA.auc(pop_bolus_sd, interval = [(0,12), (12,24)], method=:linuplogdown, sigdigits=3)


# If we want to look at a parameter for 1 individual 
thalf_4     = NCA.thalf(pop_bolus_sd[4], sigdigits=3) # Half-life calculation for 4th individual




# Run Annotated NCA for Final Report 
nca_bolus_sd_report = run_nca(pop_bolus_sd, sigdigits=3,
                        studyid="STUDY-001",
                        studytitle="Phase 1 Drug Trial: IV Bolus Single Dose", # required
                        author = [("Author 1", "Author 2")], # required
                        sponsor = "PumasAI",
                        date=Dates.now(),
                        conclabel="Drug Concentration (mg/L)",
                        timelabel="Time (hr)",
                        versionnumber=v"0.1",)
#


# Summarize Results of Interest for Final Report
param_summary_bolus_sd  = summarize(nca_bolus_sd_report.reportdf, 
                            parameters = [:half_life, 
                                          :tmax, 
                                          :cmax, 
                                          :auclast, 
                                          :vz_obs, 
                                          :cl_obs, 
                                          :aucinf_obs])


# Generate NCA Report 
report(nca_bolus_sd_report, param_summary_bolus_sd)



# Look at Individual Fits 
individual_fits = subject_fits(nca_bolus_sd,
             axis = (xlabel = "Time (hr)", 
                     ylabel = "Drug Concentration (mg/L)",
                     yscale = log10),
             separate = true, paginate = true,
             limit = 16, columns = 4, rows = 4, 
             facet = (combinelabels = true,))
#
individual_fits[1]














###################################################################
#                                                                 #
#               IV INFUSION SINGLE DOSE EXAMPLE                   #
#                                                                 #
###################################################################



# Load Data
df_inf_sd = CSV.read("WCOP//nca_data//iv_infusion_sd.csv", DataFrame, missingstrings=["NA", ".", ""])



# Map Dataframe to NCA Population
pop_inf_sd = read_nca(df_inf_sd,
                id            = :id,
                time          = :time,
                observations  = :conc,
                amt           = :amt,
                route         = :route,
                group         = [:group,],
                llq           = 0.001)
#



# Preview Data 
## individual plot - linear scale 
obsvstimes = observations_vs_time(pop_inf_sd[1])




# Perform Simple NCA
nca_inf_sd = run_nca(pop_inf_sd, sigdigits=3)






# Generate Individual NCA Parameters 
vz        = NCA.vz(pop_inf_sd, sigdigits=3)  # Volume of Distribution/F, in this case since the drug is given orally
cl        = NCA.cl(pop_inf_sd, sigdigits=3)  # Clearance/F, in this case since the drug is given orally
lambdaz   = NCA.lambdaz(pop_inf_sd, threshold=3, sigdigits=3)  # Terminal Elimination Rate Constant, threshold=3 specifies the max no. of time point used for calculation
thalf     = NCA.thalf(pop_inf_sd, sigdigits=3) # Half-life calculation for 4th individual
cmax_d    = NCA.cmax(pop_inf_sd, normalize=true, sigdigits=3) # Dose Normalized Cmax
mrt       = NCA.mrt(pop_inf_sd, sigdigits=3) # Mean residence time
aumc      = NCA.aumc(pop_inf_sd, method=:linlog, sigdigits=3) # AUMC calculation, using :linlog method
individual_params    = innerjoin(vz,cl,lambdaz,thalf,cmax_d,mrt,aumc, on=[:id,:group], makeunique=true) # include group to innerjoin**




# Other AUC calculation options 
auc0_12   = NCA.auc(pop_inf_sd, interval=(0,12), method=:linuplogdown, sigdigits=3) #various other methods are :linear, :linlog
auc12_24  = NCA.auc(pop_inf_sd, interval=(12,24), method=:linuplogdown, sigdigits=3) #looking at auc 12 to 24 hours (can make this interval anything!)
partial_aucs = NCA.auc(pop_inf_sd, interval = [(0,12), (12,24)], method=:linuplogdown, sigdigits=3)



# Run Annotated NCA for Final Report 
nca_inf_sd_report = run_nca(pop_inf_sd, sigdigits=3,
                        studyid="STUDY-002",
                        studytitle="Phase 1 Drug Trial: IV Infusion Single Dose", # required
                        author = [("Author 1", "Author 2")], # required
                        sponsor = "PumasAI",
                        date=Dates.now(),
                        conclabel="Drug Concentration (mg/L)",
                        timelabel="Time (hr)",
                        versionnumber=v"0.1",)
#


# Summarize Results of Interest for Final Report
param_summary_inf_sd  = summarize(nca_inf_sd_report.reportdf, 
                                stratify_by=[:group,], # stratifying by group so we can compare each dose 
                                parameters = [:half_life, 
                                          :tmax, 
                                          :cmax, 
                                          :auclast, 
                                          :vz_obs, 
                                          :cl_obs, 
                                          :aucinf_obs])


# Generate NCA Report 
report(nca_inf_sd_report, param_summary_inf_sd)





















###################################################################
#                                                                 #
#   ORAL MULTIPLE DOSE EXAMPLE: FIRST DOSE & STEADY STATE DATA    #
#                                                                 #
###################################################################



# Load Data
df_oral_first_ss = CSV.read("WCOP//nca_data//oral_md_first_ss.csv", DataFrame, missingstrings=["NA", ".", ""])




# Map Dataframe to NCA Population
pop_oral_first_ss = read_nca(df_oral_first_ss,
                id            = :id,
                time          = :time,
                observations  = :conc,
                amt           = :amt,
                route         = :route,
                ii            = :ii,
                ss            = :ss, 
                llq           = 0.001)
#




# Perform Simple NCA
nca_oral_first_ss = run_nca(pop_oral_first_ss, sigdigits=3)



# Generate Individual NCA Parameters 
vz        = NCA.vz(pop_oral_first_ss, sigdigits=3)  # Volume of Distribution/F, in this case since the drug is given orally
cl        = NCA.cl(pop_oral_first_ss, sigdigits=3)  # Clearance/F, in this case since the drug is given orally
lambdaz   = NCA.lambdaz(pop_oral_first_ss, threshold=3, sigdigits=3)  # Terminal Elimination Rate Constant, threshold=3 specifies the max no. of time point used for calculation
thalf     = NCA.thalf(pop_oral_first_ss, sigdigits=3) # Half-life calculation 
cmax_d    = NCA.cmax(pop_oral_first_ss, normalize=true, sigdigits=3) # Dose Normalized Cmax
mrt       = NCA.mrt(pop_oral_first_ss, sigdigits=3) # Mean residence time
aumc      = NCA.aumc(pop_oral_first_ss, method=:linlog, sigdigits=3) # AUMC calculation, using :linlog method
individual_params      = innerjoin(vz,cl,lambdaz,thalf,cmax_d,mrt,aumc, on=[:id], makeunique=true)


auc0_12   = NCA.auc(pop_inf_sd, interval=(0,12), method=:linuplogdown, sigdigits=3) #various other methods are :linear, :linlog
auc12_24  = NCA.auc(pop_inf_sd, interval=(12,24), method=:linuplogdown, sigdigits=3) #looking at auc 12 to 24 hours (can make this interval anything!)
partial_aucs = NCA.auc(pop_inf_sd, interval = [(0,12), (12,24)], method=:linuplogdown, sigdigits=3)





# Run Annotated NCA for Final Report 
nca_oral_first_ss_report = run_nca(pop_oral_first_ss, sigdigits=3,
                        studyid="STUDY-003",
                        studytitle="Phase 1 Drug Trial: Oral Multiple Dosing (First & SS)", # required
                        author = [("Author 1", "Author 2")], # required
                        sponsor = "PumasAI",
                        date=Dates.now(),
                        conclabel="Drug Concentration (mg/L)",
                        timelabel="Time (hr)",
                        versionnumber=v"0.1",)
#


# Summarize Results of Interest for Final Report
param_summary_oral_first_ss  = summarize(nca_oral_first_ss_report.reportdf, 
                                        parameters = [:half_life, 
                                                        :tmax, 
                                                        :cmax, 
                                                        :auclast, 
                                                        :vz_f_obs, # note the change because of oral 
                                                        :cl_f_obs])
#


# Generate NCA Report 
report(nca_oral_first_ss_report, param_summary_oral_first_ss)
























###################################################################
#                                                                 #
#      ORAL MULTIPLE DOSE EXAMPLE: STEADY STATE DATA ONLY         #
#                                                                 #
###################################################################



# Load Data
df_oral_ss = CSV.read("WCOP//nca_data//oral_md_ss_only.csv", DataFrame, missingstrings=["NA", ".", ""])




# Map Dataframe to NCA Population
pop_oral_ss = read_nca(df_oral_ss,
                id            = :id,
                time          = :tad,
                observations  = :conc,
                amt           = :amt,
                route         = :route,
                ii            = :ii,
                ss            = :ss, 
                llq           = 0.001)
#




# Perform Simple NCA
nca_oral_first_ss = run_nca(pop_oral_first_ss, sigdigits=3)



# Generate Individual NCA Parameters 
vz        = NCA.vz(pop_oral_first_ss, sigdigits=3)  # Volume of Distribution/F, in this case since the drug is given orally
cl        = NCA.cl(pop_oral_first_ss, sigdigits=3)  # Clearance/F, in this case since the drug is given orally
lambdaz   = NCA.lambdaz(pop_oral_first_ss, threshold=3, sigdigits=3)  # Terminal Elimination Rate Constant, threshold=3 specifies the max no. of time point used for calculation
thalf     = NCA.thalf(pop_oral_first_ss, sigdigits=3) # Half-life calculation 
cmax_d    = NCA.cmax(pop_oral_first_ss, normalize=true, sigdigits=3) # Dose Normalized Cmax
mrt       = NCA.mrt(pop_oral_first_ss, sigdigits=3) # Mean residence time
aumc      = NCA.aumc(pop_oral_first_ss, method=:linlog, sigdigits=3) # AUMC calculation, using :linlog method
individual_params      = innerjoin(vz,cl,lambdaz,thalf,cmax_d,mrt,aumc, on=[:id], makeunique=true)


auc0_12   = NCA.auc(pop_inf_sd, interval=(0,12), method=:linuplogdown, sigdigits=3) #various other methods are :linear, :linlog
auc12_24  = NCA.auc(pop_inf_sd, interval=(12,24), method=:linuplogdown, sigdigits=3) #looking at auc 12 to 24 hours (can make this interval anything!)
partial_aucs = NCA.auc(pop_inf_sd, interval = [(0,12), (12,24)], method=:linuplogdown, sigdigits=3)





# Run Annotated NCA for Final Report 
nca_oral_ss_report = run_nca(pop_oral_ss, sigdigits=3,
                        studyid="STUDY-004",
                        studytitle="Phase 1 Drug Trial: Oral Multiple Dosing (SS Only)", # required
                        author = [("Author 1", "Author 2")], # required
                        sponsor = "PumasAI",
                        date=Dates.now(),
                        conclabel="Drug Concentration (mg/L)",
                        timelabel="Time (hr)",
                        versionnumber=v"0.1",)
#


# Summarize Results of Interest for Final Report
param_summary_oral_ss  = summarize(nca_oral_ss_report.reportdf, 
                                        parameters = [:half_life, 
                                                        :tmax, 
                                                        :cmax, 
                                                        :auclast, 
                                                        :vz_f_obs, 
                                                        :cl_f_obs])
#


# Generate NCA Report 
report(nca_oral_ss_report, param_summary_oral_ss)

