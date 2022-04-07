using CSV, Pumas, PumasUtilities, Bioequivalence

# Read Data
# Mind relative paths when importing data
pkdata = CSV.read("nca_data/data4NCA_sad.csv", DataFrame)

# map NCA data from CSV
ncapop = read_nca(pkdata;
                  observations=:concentration,
                  amt=:amt,
                  route=:route,
                  group=[:doselevel])

# run the complete NCA analysis
pk_nca = run_nca(ncapop; sigdigits=2)

strata = [:doselevel]
params = [:cmax, :aucinf_obs]
output = summarize(pk_nca.reportdf;
                   stratify_by=strata,
                   parameters=params)

# Write output to CSV file
CSV.write("nca_output.csv", output)
