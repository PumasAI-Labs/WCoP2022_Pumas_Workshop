using CSV, Pumas, PumasUtilities, Bioequivalence

# Read Data from a csv file into a DataFrame
pkdata = CSV.read("be_data/data4be.csv", DataFrame)

# BE with Cmax
cmax_output = pumas_be(pkdata;
                       endpoint=:Cmax,
                       id=:id,
                       sequence=:sequence,
                       period=:period)

# BE with AUC
auc_output = pumas_be(pkdata;
                      endpoint=:AUC,
                      id=:id,
                      sequence=:sequence,
                      period=:period)

# For more information on BE in Pumas visit:
# https://docs.pumas.ai/stable/bioequivalence/be/
