# Data collection code
using XLSX, CSV, MAT
using LinearAlgebra
using DataFrames, Dates, Statistics
using Plots
include("DataCooking.jl") # Include to import some functions

# Positive data available 
raw_dfΔD = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAnewpos.csv", DataFrame)
dfΔD = raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :] # "New Detected" individuals from https://covid19.infn.it/iss/
dfΔD_2 = sum_columns(raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :])
dfΔD_tot = sum_total(raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :])
df_perc_D = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITApositivi.csv", DataFrame)

# Hospitalised data available 
raw_dfΔD = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAnewpos.csv", DataFrame)
dfΔD = raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :] # "New Detected" individuals from https://covid19.infn.it/iss/
dfΔD_2 = sum_columns(raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :])
dfΔD_tot = sum_total(raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :])
df_perc_D = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITApositivi.csv", DataFrame)
