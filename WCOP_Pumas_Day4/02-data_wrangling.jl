# Don't forget to check: https://tutorials.pumas.ai/html/DataWranglingInJulia/04-read_data.html

using CSV
using DataFramesMeta

##########################
#     I/O CSV Files      #
##########################

df = CSV.read("data/iv_sd_demogs.csv", DataFrame)

# Different delimiters and decimals
# using the keyword arguments `delim` and `decimal`
df_eu = CSV.read("data/iv_sd_demogs_eu.csv", DataFrame; delim=';', decimal=',')

# Custom types for columns
df_custom_types = CSV.read("data/iv_sd_demogs.csv", DataFrame;
	types=Dict(:ID => String,
			   :ISMALE => Bool)
)

# Selecting and dropping columns
df_select_names = CSV.read("data/iv_sd_demogs.csv", DataFrame; select=["ID", "AGE"])
df_select_idxs = CSV.read("data/iv_sd_demogs.csv", DataFrame; select=[1, 2])

df_drop_names = CSV.read("data/iv_sd_demogs.csv", DataFrame; drop=["SCR", "eGFR"])
df_drop_idx = CSV.read("data/iv_sd_demogs.csv", DataFrame; drop=[4, 6])

# Missing values
df_missing = CSV.read("data/iv_sd_demogs_missing.csv", DataFrame;     # take a look at the first row
                      missingstrings=["NA", "I don't know prof", "?"])

##########################
#     DataFramesMeta     #
##########################

# Don't forget to check: https://tutorials.pumas.ai/html/DataWranglingInJulia/05-mutating-dfmeta.html

############################################################################################
# dplyr has 50+ functions `mutate`, `mutate_if`, `mutate_at`, `mutate_all`, `rename_with`, #
# `transmute`, `transmute_if`, `transmute_at`, `transmute_all`                             #
############################################################################################

##############################################
# DataFramesMeta has 6 macros and 1 function #
##############################################

## @[r]select[!] and @[r]transform[!]
@select df :ID :AGE

# We can also use `Not()`, `Between()`, RegEx
@select df $(Between(:SCR, :eGFR))
@select df $(Not(:SCR, :eGFR))
@select df $(r"R$")  # ending with `R`

@select df begin
	:ID
	:eGFR_z = begin
		μ = mean(:eGFR)
		σ = std(:eGFR)
		[(x - μ) / σ for x in :eGFR]
	end
end

@transform df begin
	:eGFR_z = begin
		μ = mean(:eGFR)
		σ = std(:eGFR)
		[(x - μ) / σ for x in :eGFR]
	end
end

@transform df :AGE_log = log(:AGE)  # this errors! Why?!
@rtransform df :AGE_log = log(:AGE) # this doesn't! Why?!


## @[r]subset[!]
@rsubset df :eGFR > 100
@rsubset df :eGFR > 100 :AGE < 30

@subset df :SCR .> mean(:SCR) # Why `@rsubset` would fail here?

@subset df begin
	:SCR .> mean(:SCR)
	:eGFR .< median(:eGFR)
end

## @orderby
@orderby df :eGFR  # ascending by default
@orderby df -:eGFR # descending

@orderby df -:ISMALE :eGFR  # several conditions

## @chain
# This is the pipe!
@chain df begin
    @select $(Between(:SCR, :eGFR))
    @transform begin
        :eGFR_z = begin
            μ = mean(:eGFR)
            σ = std(:eGFR)
            [(x - μ) / σ for x in :eGFR]
        end
    end
    @rsubset :eGFR > 100
    @orderby -:ISMALE :eGFR
end

## groupby and @combine
@chain df begin
    groupby(:ISMALE)
    @combine begin
        :AGE_μ = mean(:AGE)
        :WEIGHT_μ = mean(:WEIGHT)
        :total = length(:ID)
        :high_eGFR = count(>(80), :eGFR)
    end
end

# using a lazy evaluation `$()` with broadcasting `.`
@chain df begin
    groupby(:ISMALE)
    @combine $([:AGE, :WEIGHT] .=> [mean median])
end
