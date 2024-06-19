# Code to plot dataset results

df_H = CSV.read("Reconstructed Datasets/Reconstructed_H.csv", DataFrame)
df_D = CSV.read("Reconstructed Datasets/Reconstructed_D.csv", DataFrame)
dfΔD_2 = CSV.read("Reconstructed Datasets/NewPositive_ΔD.csv", DataFrame)  
dfΔE_2 = CSV.read("Reconstructed Datasets/NewDeceased_ΔE.csv", DataFrame)  
DailyTrend_df = CSV.read("/Users/marcodelloro/Desktop/Pandemic-System-Identification/National Trends/DailyTrend_ITA.csv", DataFrame)
E_dy = DailyTrend_df.deceduti

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=true) 
scalefontsizes(1)

# ΔI values compared to ΔD (age compartment based)
plot(I_age.u40.Date, I_age.u40.ΔI, ribbon=I_age.u40.varΔI,fillalpha=.5,label="ΔI-Capture Recapture", xlabel="Date", ylabel="Population", title=" ΔI - under 40", linewidth=1.5)
scatter!(dfΔD_2.Date, dfΔD_2.u40, label="ΔD-Data", ms=2, mc=:red, markerstrokecolor=:red, ma=0.6)

# I values compared to D (age compartment based)
plot(I_age.u40.Date, I_age.u40.I, label="I-Capture Recapture", xlabel="Date", ylabel="Population", title=" I - under 40", linewidth=2)
scatter!(dfΔD_2.Date, df_D.u40, label="D-Data", ms=2, mc=:red, markerstrokecolor=:red, ma=0.6)

plot(I_age.mid.Date, I_age.mid.I,label="I-Capture Recapture", xlabel="Date", ylabel="Population", title=" I - 40÷60", linewidth=2)
scatter!(dfΔD_2.Date, df_D.mid, label="D-Data", ms=2, mc=:red, markerstrokecolor=:red, ma=0.6)

plot(I_age.mid.Date, I_age.old.I,label="I-Capture Recapture", xlabel="Date", ylabel="Population", title=" I - 60÷80", linewidth=2)
scatter!(dfΔD_2.Date, df_D.old, label="D-Data", ms=2, mc=:red, markerstrokecolor=:red, ma=0.6)

plot(I_age.mid.Date, I_age.ger.I,label="I-Capture Recapture", xlabel="Date", ylabel="Population", title=" I - 80+", linewidth=2)
scatter!(dfΔD_2.Date, df_D.ger, label="D-Data", ms=2, mc=:red, markerstrokecolor=:red, ma=0.6)

# H values between sources
plot(df_H.Date, df_H_tot.Tot, label="AgeGroups Sum", linewidth=2.5)

# H values for the different age groups vs SUM
plot(df_H.Date, df_H.u40, label="under 40", xlabel="Date", ylabel="Population", title=" H - Healed Individuals", linewidth=1.5)
plot!(df_H.Date, df_H.mid, label="40÷60", linewidth=1.5)
plot!(df_H.Date, df_H.old, label="60÷80", linewidth=1.5)
plot!(df_H.Date, df_H.ger, label="80+", linewidth=1.5)
plot!(df_H.Date, df_H_tot.Tot, label="Age Groups Sum", linewidth=2.5)
