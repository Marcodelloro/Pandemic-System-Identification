# Data collection code
using XLSX, CSV, MAT
using LinearAlgebra
using DataFrames, Dates, Statistics
using Plots

start_date = Date("2020-08-31")
end_date = Date("2021-10-03")
df_H = CSV.read("Reconstructed_H.csv", DataFrame)

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

#                                   ---- Data divided by age based on https://covid19.infn.it/iss/ Database ---- 
# Positive d ata available 
raw_dfΔD = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age New Cases Data/issITAnewpos.csv", DataFrame)
dfΔD = raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :] 
dfΔD_2 = sum_columns(raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :])
dfΔD_tot = sum_total(raw_dfΔD[(raw_dfΔD.data .>= start_date) .& (raw_dfΔD.data .<= end_date), :])
df_perc_D = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITApositivi.csv", DataFrame)

# Hospitalised data available 
raw_dfΔT1 = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age New Cases Data/issITAnew_ricoveri.csv", DataFrame)
dfΔT1 = raw_dfΔT1[(raw_dfΔT1.data .>= start_date) .& (raw_dfΔT1.data .<= end_date), :] 
dfΔDT1_2 = sum_columns(raw_dfΔT1[(raw_dfΔT1.data .>= start_date) .& (raw_dfΔT1.data .<= end_date), :])
dfΔT1_tot = sum_total(raw_dfΔT1[(raw_dfΔT1.data .>= start_date) .& (raw_dfΔT1.data .<= end_date), :])
df_perc_T1 = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAhosp.csv", DataFrame)

# Intensive Care data available 
raw_dfΔT2 = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age New Cases Data/issITAnew_ICU.csv", DataFrame)
dfΔT2 = raw_dfΔT2[(raw_dfΔT2.data .>= start_date) .& (raw_dfΔT2.data .<= end_date), :]
dfΔT2_2 = sum_columns(raw_dfΔT2[(raw_dfΔT2.data .>= start_date) .& (raw_dfΔT2.data .<= end_date), :])
dfΔT2_tot = sum_total(raw_dfΔT2[(raw_dfΔT2.data .>= start_date) .& (raw_dfΔT2.data .<= end_date), :])
df_perc_T2 = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAicu.csv", DataFrame)

# Deceased/Deaths data available 
raw_dfΔE = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age New Cases Data/issITAnewdeceased.csv", DataFrame)
dfΔE = sum_columns(raw_dfΔE)
dfΔE_tot = sum_total(raw_dfΔE)
dfΔE_tot[!,"Cumulated_Sum"]= cumsum(dfΔE_tot.Tot) # total sum based on ISS data
dfΔE_2 = dfΔE_tot[(dfΔE_tot.Date .>= start_date) .& (dfΔE_tot.Date .<= end_date), :]
df_perc_E = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/Age Compatments Distributions/issITAdeceased.csv", DataFrame)

# Healed Individuals data available 
df_H_tot = sum_total(df_H)

#                                 ---- National Overall Data based on https://github.com/pcm-dpc/COVID-19 ---- 

# Positive data available 
DailyTrend_df = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/National Trends/DailyTrend_ITA.csv", DataFrame)
WeekTrend_df = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/National Trends/WeekTrend_ITA.csv", DataFrame)
DailyTrend_df = DailyTrend_df[(DailyTrend_df.data .>= start_date) .& (DailyTrend_df.data .<= end_date), :]
WeekTrend_df = WeekTrend_df[(WeekTrend_df.giorno .>= start_date) .& (WeekTrend_df.giorno .<= end_date), :]

D_wk = WeekTrend_df.mm_totale_positivi # Currently positive individuals
E_wk = WeekTrend_df.mm_deceduti # Currently deceased individuals
H_wk = WeekTrend_df.mm_guariti # Currently Healed individuals
Dtot_wk =  WeekTrend_df.mm_totale_casi # Cumulated sum of cases up to date 

D_dy = DailyTrend_df.totale_positivi # Currently positive individuals
ΔD_dy = DailyTrend_df.nuovi_positivi # New positive individuals (overall)
E_dy = DailyTrend_df.deceduti # Currently Deceased Individuals
ΔE_dy =  DailyTrend_df.diff_deceduti # New deceased individuals (overall)

T2_dy = DailyTrend_df.terapia_intensiva # Current ICUs numbers
T1_dy = DailyTrend_df.ricoverati # Current Hospitalised individuals 
ΔT1_dy = DailyTrend_df.diff_ricoverati # Difference (+;-) of ICUs in consectuive days
ΔT2_dy = DailyTrend_df.diff_terapia_intensiva # Current Hospitalised individuals 

H_dy =  DailyTrend_df.guariti # Currently Healed individuals 
ΔH_dy = DailyTrend_df.guariti # New daily Healed Individuals

#                                 ---- Comparison plots ---- 

# D values for the different age groups vs SUM
plot(dfΔD_2.Date, dfΔD_2.u40, label="u40", xlabel="Date", ylabel="Population", title=" D - Detected Individuals", linewidth=1.5)
plot!(dfΔD_2.Date, dfΔD_2.mid, label="mid", linewidth=1.5)
plot!(dfΔD_2.Date, dfΔD_2.old, label="old", linewidth=1.5)
plot!(dfΔD_2.Date, dfΔD_2.ger, label="ger", linewidth=1.5)
plot!(dfΔD_2.Date, dfΔD_tot.Tot, label="AgeGroups Sum", linewidth=2.5)

# D values between sources
plot(dfΔD_2.Date, dfΔD_tot.Tot, label="Tot Age Groups Sum", xlabel="Date", ylabel="Population", title=" D - Detected Individuals", linewidth=1.5)
plot!(dfΔD_2.Date, D_wk, label="Weekly Mean Currently D", linewidth=1.5)
plot!(dfΔD_2.Date, Dtot_wk, label="Weekly Cumulated", linewidth=1.5)
#plot!(dfΔD_2.Date, D_dy, label="Daily Currently D", linewidth=1.5)
plot!(dfΔD_2.Date, ΔD_dy, label="New Daily", linewidth=1.5)

# E values between sources
plot(dfΔD_2.Date, dfΔE_2.Cumulated_Sum, label="Tot Age Groups Sum", xlabel="Date", ylabel="Population", title="E - Deceased", linewidth=2.5)
plot!(dfΔD_2.Date, E_dy, label="Daily Values", xlabel="Date", ylabel="Population", title="E - Deceased", linewidth=2.5)

# H values for the different age groups vs SUM
plot(df_H.Date, df_H.u40, label="u40", xlabel="Date", ylabel="Population", title=" D - Detected Individuals", linewidth=1.5)
plot!(df_H.Date, df_H.mid, label="mid", linewidth=1.5)
plot!(df_H.Date, df_H.old, label="old", linewidth=1.5)
plot!(df_H.Date, df_H.ger, label="ger", linewidth=1.5)
plot!(df_H.Date, df_H_tot.Tot, label="AgeGroups Sum", linewidth=2.5)

# H values between sources
plot(df_H.Date, df_H_tot.Tot, label="AgeGroups Sum", linewidth=2.5)
plot!(df_H.Date, H_dy, label="Daily Values", xlabel="Date", ylabel="Population", title="H - Healed Individuals", linewidth=2.5)
