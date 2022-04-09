# Don't forget to check: https://tutorials.pumas.ai/html/DataWranglingInJulia/04-read_data.html

using CSV
using DataFramesMeta

##########################
#     I/O CSV Files      #
##########################

df = CSV.read("data/iv_sd_demogs.csv", DataFrame)

1. Selecting and dropping columns
   2. Different delimiters and decimals