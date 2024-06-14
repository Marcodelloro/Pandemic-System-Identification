using XLSX, CSV, MAT
using LinearAlgebra
using DataFrames, Dates, Statistics
using Plots

# Age percentage in compartments from iss database loading
dfD = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITApositivi.csv", DataFrame)
dfT1 = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAhosp.csv", DataFrame)
dfT2 = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAicu.csv", DataFrame)
dfE = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAdeceased.csv", DataFrame)

start_date = Date("2020-08-31")
end_date = Date("2021-10-03")
date = collect(start_date:end_date)
N_yng = 23205875 # Total under40 pop
N_mid = 18142711 # Total 40-60 pop
N_old = 13408810 # Total 60-80 pop
N_eld = 3951057  # Total 80+ pop

function df_normrows(df::DataFrame)
    for i in 1:nrow(df)
        row_sum = sum(skipmissing(df[i, 2:end]))
        for j in 2:ncol(df)
           df[i, j] = df[i, j] / row_sum
        end
    end
    return df
end

dfD_filt = df_normrows(dfD[(dfD.data .>= start_date) .& (dfD.data .<= end_date), :])
dfT1_filt = df_normrows(dfT1[(dfT1.data .>= start_date) .& (dfT1.data .<= end_date), :])
dfT2_filt = df_normrows(dfT2[(dfT2.data .>= start_date) .& (dfT2.data .<= end_date), :])
dfE_filt = df_normrows(dfE[(dfE.data .>= start_date) .& (dfE.data .<= end_date), :])

# Contact matrix in different locations
c_ovr = XLSX.readdata("ContactMarticesRaw.xlsx", "Overall!A1:P16")
c_home = XLSX.readdata("ContactMarticesRaw.xlsx", "Home!A1:P16")
c_sch= XLSX.readdata("ContactMarticesRaw.xlsx", "School!A1:P16")
c_work = XLSX.readdata("ContactMarticesRaw.xlsx", "Work!A1:P16")
c_othr = XLSX.readdata("ContactMarticesRaw.xlsx", "Others!A1:P16")

# 3x3 matrix rearranging to match available data
x = [ ones(1,8) zeros(1,8); 
      zeros(1,8) ones(1,4) zeros(1,4);
      zeros(1,8) zeros(1,4) ones(1,4)   ]
c_ovr = x*c_ovr*x' 
c_home = x*c_home*x'
c_sch= x*c_sch*x'
c_work = x*c_work*x'
c_othr =x*c_othr*x'

# Loading of the Italian data from SIDTHE
SIDTTHEfile = matopen("SIDTTHE_data_DEF.mat")
SIDTTHE_data = read(SIDTTHEfile, "SIDTTHE_data")
close(SIDTTHEfile)

# Struct with the data for the different compartments
struct Compartment
    name::String
    young::Matrix{Float64}
    mid::Matrix{Float64}
    old::Matrix{Float64}
    ger::Matrix{Float64}
end

D_data = Compartment("Detected", SIDTTHE_data[2, 1]["data"] .* dfD_filt[!,"0-39 anni"]', SIDTTHE_data[2, 1]["data"] .* dfD_filt[!,"40-59 anni"]',SIDTTHE_data[2, 1]["data"] .* dfD_filt[!,"60-79 anni"]', (SIDTTHE_data[2, 1]["data"] .* dfD_filt[!,"80-89 anni"]') + (SIDTTHE_data[2, 1]["data"] .* dfD_filt[!,"≥90 anni"]'))
T1_data = Compartment("Hospitalised", SIDTTHE_data[3, 1]["data"] .* dfT1_filt[!,"0-39 anni"]', SIDTTHE_data[3, 1]["data"] .* dfT1_filt[!,"40-59 anni"]',SIDTTHE_data[3, 1]["data"] .* dfT1_filt[!,"60-79 anni"]', (SIDTTHE_data[3, 1]["data"] .* dfT1_filt[!,"80-89 anni"]') + (SIDTTHE_data[3, 1]["data"] .* dfT1_filt[!,"≥90 anni"]'))
T2_data = Compartment("IntensiveCare", SIDTTHE_data[4, 1]["data"] .* dfT2_filt[!,"0-39 anni"]', SIDTTHE_data[4, 1]["data"] .* dfT2_filt[!,"40-59 anni"]',SIDTTHE_data[4, 1]["data"] .* dfT2_filt[!,"60-79 anni"]', (SIDTTHE_data[4, 1]["data"] .* dfT2_filt[!,"80-89 anni"]') + (SIDTTHE_data[4, 1]["data"] .* dfT2_filt[!,"≥90 anni"]'))
E_data = Compartment("Deceased", SIDTTHE_data[6, 1]["data"] .* dfE_filt[!,"0-39 anni"]', SIDTTHE_data[6, 1]["data"] .* dfE_filt[!,"40-59 anni"]',SIDTTHE_data[6, 1]["data"] .* dfE_filt[!,"60-79 anni"]', (SIDTTHE_data[6, 1]["data"] .* dfE_filt[!,"80-89 anni"]') + (SIDTTHE_data[6, 1]["data"] .* dfE_filt[!,"≥90 anni"]'))

D_df = DataFrame( Date = date, u40 = vec(D_data.young), mid = vec(D_data.mid),
    old = vec(D_data.old), ger = vec(D_data.ger)
)

T1_df = DataFrame( Date = date, u40 = vec(T1_data.young), mid = vec(T1_data.mid),
    old = vec(T1_data.old), ger = vec(T1_data.ger)
)

T2_df = DataFrame( Date = date, u40 = vec(T2_data.young), mid = vec(T2_data.mid),
    old = vec(T2_data.old), ger = vec(T2_data.ger)
)

E_df = DataFrame( Date = date, u40 = vec(E_data.young), mid = vec(E_data.mid),
    old = vec(E_data.old), ger = vec(E_data.ger)
)

# Estimation of known Healed individuals
raw_dfΔD = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age New Cases Data/issITAnewpos.csv", DataFrame)
raw_dfΔE = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age New Cases Data/issITAnewdeceased.csv", DataFrame)

# Function to sum elements in different columns for the "New Detected" data
function sum_columns(df::DataFrame)
    result_df = DataFrame(
        Date = df[!, 1],
        u40 = Vector{Float64}(undef, nrow(df)),
        mid = Vector{Float64}(undef, nrow(df)),
        old = Vector{Float64}(undef, nrow(df)),
        ger = Vector{Float64}(undef, nrow(df))
    )
    for i in 1:nrow(df)
        result_df.u40[i] = sum(skipmissing(df[i, 2:6]); init=0)   # u40 summed cols
        result_df.mid[i] = sum(skipmissing(df[i, 7:8]); init=0)   # 40-60 summed cols
        result_df.old[i] = sum(skipmissing(df[i, 9:10]); init=0)  # 60-80 summed cols
        result_df.ger[i] = sum(skipmissing(df[i, 11:12]); init=0) # 80+ summed cols
    end

    return result_df
end

# Function to sum all the age groups in a single column
function sum_total(df::DataFrame)
    result_df = DataFrame(
        Date = df[!, 1],
        Tot = Vector{Float64}(undef, nrow(df))
    )
    for i in 1:nrow(df)
        result_df.Tot[i] = sum(skipmissing(df[i, 2:end]); init=0) # 80+ summed cols
    end

    return result_df
end

dfΔD = sum_columns(raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :])

cumsum_ΔD = DataFrame(
    Date = raw_dfΔD.data,
    u40 = cumsum(sum_columns(raw_dfΔD).u40),
    mid = cumsum(sum_columns(raw_dfΔD).mid),
    old = cumsum(sum_columns(raw_dfΔD).old),
    ger = cumsum(sum_columns(raw_dfΔD).ger)
)

Ẽ = DataFrame(
Date = raw_dfΔE.data,
    u40 = cumsum(sum_columns(raw_dfΔE).u40),
    mid = cumsum(sum_columns(raw_dfΔE).mid),
    old = cumsum(sum_columns(raw_dfΔE).old),
    ger = cumsum(sum_columns(raw_dfΔE).ger)
)

Ẽ_df = Ẽ[(Ẽ.Date .>= start_date) .& (Ẽ.Date .<= end_date), :] # New actually utilised Deceased Number (More reliable)

function reconH(df::DataFrame, csum_df::DataFrame)
    result_H = DataFrame(
        Date = df[!, 1],
        u40 = Vector{Float64}(undef, nrow(df)),
        mid = Vector{Float64}(undef, nrow(df)),
        old = Vector{Float64}(undef, nrow(df)),
        ger = Vector{Float64}(undef, nrow(df))
    )
    for j in 1:nrow(dfΔD)
        idx = findfirst(row -> row.Date == start_date, eachrow(csum_df))-1
        result_H.u40[j] = csum_df.u40[idx + j] - D_df.u40[j] - T1_df.u40[j] - T2_df.u40[j] - Ẽ_df.u40[j]
        result_H.mid[j] = csum_df.mid[idx + j] - D_df.mid[j] - T1_df.mid[j] - T2_df.mid[j] - Ẽ_df.mid[j]
        result_H.old[j] = csum_df.old[idx + j] - D_df.old[j] - T1_df.old[j] - T2_df.old[j] - Ẽ_df.old[j]
        result_H.ger[j] = csum_df.ger[idx + j] - D_df.ger[j] -  T1_df.ger[j] - T2_df.ger[j] - Ẽ_df.ger[j]
    end
    return result_H
end

actualH = reconH(dfΔD,cumsum_ΔD)
CSV.write("Reconstructed_H.csv", actualH)

