using MAT, CSV
using DataFrames
using Dates
using OrdinaryDiffEq, LinearAlgebra
using Optimization, OptimizationPolyalgorithms, OptimizationOptimJL
using Zygote, ForwardDiff
using Plots, LaTeXStrings
using Statistics, StatsBase 
using OptimizationMOI, Ipopt

S_data = CSV.read("Reconstructed Datasets/Reconstructed_S.csv", DataFrame)
I_raw =  CSV.read("Reconstructed Datasets/Infected_I.csv", DataFrame)
D_data = CSV.read("Reconstructed Datasets/Reconstructed_D.csv", DataFrame)
T1_data = CSV.read("Reconstructed Datasets/Reconstructed_T1.csv", DataFrame)
T2_data = CSV.read("Reconstructed Datasets/Reconstructed_T2.csv", DataFrame)
H_data = CSV.read("Reconstructed Datasets/Reconstructed_H.csv", DataFrame)
E_data = CSV.read("Reconstructed Datasets/Reconstructed_E.csv", DataFrame)
N_u40 = 23205875; N_mid = 18142711; 
N_old = 13408810; N_ger = 3951057;
N_mhe = 21;  N = 399; Npop = N_u40 + N_mid + N_old + N_ger; 
N_groups= [N_u40, N_mid, N_old, N_ger]

# Manipulation of raw df 
#=

I_data = vcat(I_raw[1:1, :], I_raw)
T_data = DataFrame( Date = T1_data.Date, u40 = T1_data.u40 + T2_data.u40, mid = T1_data.mid + T2_data.mid,
    old = T1_data.old + T2_data.old, ger = T1_data.ger + T2_data.ger)

ymeas_u40 = DataFrame(S_u40 = S_data.u40, I_u40 = vec(I_data.u40/Npop), D_u40 = vec(D_data.u40/Npop), T_u40 = vec(T_data.u40/Npop), H_u40 = vec(H_data.u40/Npop), E_u40 = vec(E_data.u40/Npop) )
ymeas_mid = DataFrame(S_mid = S_data.mid, I_mid = vec(I_data.mid/Npop), D_mid = vec(D_data.mid/Npop), T_mid = vec(T_data.mid/Npop), H_mid = vec(H_data.mid/Npop), E_mid = vec(E_data.mid/Npop) )
ymeas_old = DataFrame(S_old = S_data.old, I_old = vec(I_data.old/Npop), D_old = vec(D_data.old/Npop), T_old = vec(T_data.old/Npop), H_old = vec(H_data.old/Npop), E_old = vec(E_data.old/Npop) )
ymeas_ger = DataFrame(S_ger = S_data.ger, I_ger = vec(I_data.ger/Npop), D_ger = vec(D_data.ger/Npop), T_ger = vec(T_data.ger/Npop), H_ger = vec(H_data.ger/Npop), E_ger = vec(E_data.ger/Npop) )

=#

# "DomaTable" data - Reconstruct the data with the one optimized in MATLAB
DomaData = CSV.read("Reconstructed Datasets/DomaTable.csv", DataFrame)
ymeas_u40 = DataFrame(S_u40 = DomaData."S-u40", I_u40 = DomaData."I-u40" , D_u40 = DomaData."D-u40", T_u40 = DomaData."T-u40", H_u40 = DomaData."H-u40", E_u40 = DomaData."E-u40" )
ymeas_mid = DataFrame(S_mid = DomaData."S-mid", I_mid = DomaData."I-mid" , D_mid = DomaData."D-mid", T_mid = DomaData."T-mid", H_mid = DomaData."H-mid", E_mid = DomaData."E-mid" )
ymeas_old = DataFrame(S_old = DomaData."S-old", I_old = DomaData."I-old" , D_old = DomaData."D-old", T_old = DomaData."T-old", H_old = DomaData."H-old", E_old = DomaData."E-old" )
ymeas_ger = DataFrame(S_ger = DomaData."S-ger", I_ger = DomaData."I-ger" , D_ger = DomaData."D-ger", T_ger = DomaData."T-ger", H_ger = DomaData."H-ger", E_ger = DomaData."E-ger" )

# Define a struct for contact matrices from POLYMOD study

mutable struct ContactMatrices
    home::Matrix{Float64}
    schl::Matrix{Float64}
    work::Matrix{Float64}
    othr::Matrix{Float64}
end

c_struct = ContactMatrices(
    [ 19.7896   7.35764  1.43232   1.3;
      4.73514   5.62234  0.445607  0.44;
      3.26361   2.32613  4.3079    4.00; 
      3.00      2.00     4.00      3.00  ],

    [ 24.9874    2.216       0.0889567   0.0389567;
      8.94216    1.1419      0.0793523   0.0293523;
      0.666502   0.0661153   0.054174    0.014174;
      0.0566502  0.0361153   0.014174    0.004174 ],

    [ 13.2266     8.39663   0.0329148    0.0229148;
      9.40362     8.72244   0.0357202    0.0157202;
      0.0930825   0.154919  0.000840544  0.00040544;
      0.030825    0.015491  0.00040544   0.00000001 ],

    [ 44.1768     12.4868   4.61979   4.61979;
      11.3875     8.96874   4.04791   4.04791;
      6.62855     9.04265   8.95355   8.00000;
      3.62855     6.04265   9.00000   9.95355 ]
)


#= C = c_home + c_schl + c_work + c_othr =#

contfields = [:home, :schl, :work, :othr] 

for field in contfields
    c_mat = getfield(c_struct, field)
    for i = 1:size(c_mat, 2)
        for j = 1:size(c_mat, 1)
            # Symmetrization of the matrices
            c_mat[j, i] = (1 / (2 * N_groups[i])) * (c_mat[j, i] * N_groups[i] + c_mat[i, j] * N_groups[j])
        end
    end
    
    setfield!(c_struct, field, c_mat)
end

c_tot = c_struct.home + c_struct.schl + c_struct.work + c_struct.othr
normcoef = sum(c_tot, dims=2) 
# Guessed values of lambda, assuming that older people heal slower
λ1 = 1/10; λ2 = 1/15; λ3= 1/21; λ4 = 1/30; 

C = c_struct.home ./ normcoef

function SIDTHEage!(du, u, p, t)
    # State variables
    s1, i1, d1, t1, h1, e1,
    s2, i2, d2, t2, h2, e2,
    s3, i3, d3, t3, h3, e3,
    s4, i4, d4, t4, h4, e4 = u 
    
    # Parameters
    α1, α2, α3, α4,
    γ1, γ2, γ3, γ4,
    δ1, δ2, δ3, δ4,
    σ1, σ2, σ3, σ4,
    τ1, τ2, τ3, τ4 = p

    # Set of equations for u40 (group 1)
    du[1]  = -s1 * α1 * (C[1, 1] * i1 + C[1, 2] * i2 + C[1, 3] * i3 + C[1, 4] * i4)                                   # dS1/dt
    du[2]  = s1 * α1 * (C[1, 1] * i1 + C[1, 2] * i2 + C[1, 3] * i3 + C[1, 4] * i4) - ((λ1*γ1)/(λ1+γ1)) * i1  - γ1*i1  # dI1/dt
    du[3]  = i1 * γ1 - d1 * (λ1 + δ1)                                                                                 # dD1/dt
    du[4]  = δ1 * d1 - (σ1 + τ1) * t1                                                                                 # dT1/dt
    du[5]  = ((λ1*γ1)/(λ1+γ1)) * i1 + d1 * λ1 + t1 * σ1                                                          # dH1/dt
    du[6]  = t1 * τ1                                                                                                  # dE1/dt

    # Set of equations for 40-60 (middle aged) (group 2)
    du[7]  = -s2 * α2 * (C[2, 1] * i1 + C[2, 2] * i2 + C[2, 3] * i3 + C[2, 4] * i4)                                   # dS2/dt
    du[8]  = s2 * α2 * (C[2, 1] * i1 + C[2, 2] * i2 + C[2, 3] * i3 + C[2, 4] * i4) - ((λ2*γ2)/(λ2+γ2)) * i2  - γ2*i2  # dI2/dt
    du[9]  = i2 * γ2 - d2 * (λ2 + δ2)                                                                                 # dD2/dt
    du[10] = δ2 * d2 - (σ2 + τ2) * t2                                                                                 # dT2/dt
    du[11] = ((λ2*γ2)/(λ2+γ2)) * i2 + d2 * λ2 + t2 * σ2                                                          # dH2/dt
    du[12] = t2 * τ2                                                                                                  # dE2/dt

    # Set of equations for 60-80 (old aged) (group 3)
    du[13] = -s3 * α3 * (C[3, 1] * i1 + C[3, 2] * i2 + C[3, 3] * i3 + C[3, 4] * i4)                                    # dS3/dt
    du[14] = s3 * α3 * (C[3, 1] * i1 + C[3, 2] * i2 + C[3, 3] * i3 + C[3, 4] * i4) - ((λ3*γ3)/(λ3+γ3)) * i3  - γ3*i3   # dI3/dt
    du[15] = i3 * γ3 - d3 * (λ3 + δ3)                                                                                  # dD3/dt
    du[16] = δ3 * d3 - (σ3 + τ3) * t3                                                                                  # dT3/dt
    du[17] = ((λ3*γ3)/(λ3+γ3)) * i3 + d3 * λ3 + t3 * σ3                                                           # dH3/dt
    du[18] = t3 * τ3                                                                                                   # dE3/dt

    # Set of equations for 80+ (group 4)
    du[19] = -s4 * α4 * (C[4, 1] * i1 + C[4, 2] * i2 + C[4, 3] * i3 + C[4, 4] * i4)                                    # dS4/dt
    du[20] = s4 * α4 * (C[4, 1] * i1 + C[4, 2] * i2 + C[4, 3] * i3 + C[4, 4] * i4) - ((λ4*γ4)/(λ4+γ4)) * i4   - γ4*i4  # dI4/dt
    du[21] = i4 * γ4 - d4 * (λ4 + δ4)                                                                                  # dD4/dt
    du[22] = δ4 * d4 - (σ4 + τ4) * t4                                                                                  # dT4/dt
    du[23] = ((λ4*γ4)/(λ4+γ4)) * i4 + d4*λ4 + t4 * σ4                                                             # dH4/dt
    du[24] = t4 * τ4                                                                                                   # dE4/dt
end



# Initial conditions for the simulation 
u0 = [ymeas_u40.S_u40[1], ymeas_u40.I_u40[1], ymeas_u40.D_u40[1], ymeas_u40.T_u40[1], ymeas_u40.H_u40[1], ymeas_u40.E_u40[1], 
      ymeas_mid.S_mid[1], ymeas_mid.I_mid[1], ymeas_mid.D_mid[1], ymeas_mid.T_mid[1], ymeas_mid.H_mid[1], ymeas_mid.E_mid[1], 
      ymeas_old.S_old[1], ymeas_old.I_old[1], ymeas_old.D_old[1], ymeas_old.T_old[1], ymeas_old.H_old[1], ymeas_old.E_old[1],  
      ymeas_ger.S_ger[1], ymeas_ger.I_ger[1], ymeas_ger.D_ger[1], ymeas_ger.T_ger[1], ymeas_ger.H_ger[1], ymeas_ger.E_ger[1] ]

# First test parameters
#= α01 = 3;           α02 = 3;          α03 = 3;       α04 = 3
γ01 = 0.1;         γ02 = 0.05;       γ03 = 0.02;    γ04 = 0.05
δ01 = 0.8e-3;      δ02 = 2.5e-3;     δ03 = 5e-3;    δ04 = 0.5
σ01 = 0.07;        σ02 = 0.06;       σ03 = 0.04;    σ04 = 0.04
τ01 = 1.4e-3;      τ02 = 5e-3;       τ03 = 0.013;   τ04 = 0.03 =#

# "DomaTable" test parameters
α01 = DomaData."alp-u40"[1];       α02 = DomaData."alp-mid"[1];          
α03 = DomaData."alp-old"[1];       α04 = DomaData."alp-ger"[1]

γ01 = DomaData."gam-u40"[1];       γ02 = DomaData."gam-mid"[1];      
γ03 = DomaData."gam-old"[1];       γ04 = DomaData."gam-ger"[1]

δ01 = DomaData."dlt-u40"[1];       δ02 = DomaData."dlt-mid"[1];  
δ03 = DomaData."dlt-old"[1];       δ04 = DomaData."dlt-ger"[1];

σ01 = DomaData."sgm-u40"[1];       σ02 = DomaData."sgm-mid"[1];
σ03 = DomaData."sgm-old"[1];       σ04 = DomaData."sgm-ger"[1];

τ01 = DomaData."tau-u40"[1];       τ02 = DomaData."tau-mid"[1];
τ03 = DomaData."tau-old"[1];       τ04 = DomaData."tau-ger"[1];


p0 = [α01, α02, α03, α04,
      γ01, γ02, γ03, γ04,
      δ01, δ02, δ03, δ04,
      σ01, σ02, σ03, σ04,
      τ01, τ02, τ03, τ04 ]

horizon = 28 # 4 weeks horizon
tspan = (1.0, horizon)
tsteps = collect(range(1, stop = horizon, length = horizon))
prob = ODEProblem(SIDTHEage!, u0, tspan, p0)
solsim = solve(prob, Tsit5(), saveat = tsteps)
reduce(hcat,solsim.u)


# I - Infected population
p_i1 = plot(T1_data.Date[1:horizon], [solsim[2,:], ymeas_u40.I_u40[1:horizon]],
     label=["Simulated Curve" "Measured Data"], title="I - u40", ylabel="Population", legend=:topright)
p_i2 = plot(T1_data.Date[1:horizon], [solsim[8,:], ymeas_mid.I_mid[1:horizon]],
     title="I - 40-60", ylabel="Population", legend =false)
p_i3 = plot(T1_data.Date[1:horizon], [solsim[14,:], ymeas_old.I_old[1:horizon]],
     title="I - 60-80", ylabel="Population", legend =false)
p_i4 = plot(T1_data.Date[1:horizon], [solsim[20,:], ymeas_ger.I_ger[1:horizon]],
     title="I - 80+", ylabel="Population", legend =false)
i_plot = plot(p_i1, p_i2, p_i3, p_i4,size=(900, 600))

I_sim_plot = plot(T1_data.Date[1:horizon], [solsim[2,:],solsim[8,:],solsim[14,:],solsim[20,:]],  label=["Under 40" "Middle-aged (40-60)" "Senior (60-80)" "Geriatric (80+)"], title="I - Infected", ylabel="Population", legend=:right)

# D - Infected population
p_d1 = plot(T1_data.Date[1:horizon], [solsim[3,:], ymeas_u40.D_u40[1:horizon]],
     label=["Simulated Curve" "Measured Data"], title="D - u40", ylabel="Population", legend=:topright)
p_d2 = plot(T1_data.Date[1:horizon], [solsim[9,:], ymeas_mid.D_mid[1:horizon]],
     title="D - 40-60", ylabel="Population", legend =false)
p_d3 = plot(T1_data.Date[1:horizon], [solsim[15,:], ymeas_old.D_old[1:horizon]],
     title="D - 60-80", ylabel="Population", legend =false)
p_d4 = plot(T1_data.Date[1:horizon], [solsim[21,:], ymeas_ger.D_ger[1:horizon]],
     title="D - 80+", ylabel="Population", legend =false)
d_plot = plot(p_d1, p_d2, p_d3, p_d4,size=(900, 600))

D_sim_plot = plot(T1_data.Date[1:horizon], [solsim[3,:],solsim[9,:],solsim[15,:],solsim[21,:]],  label=["Under 40" "Middle-aged (40-60)" "Senior (60-80)" "Geriatric (80+)"], title="D - Detected", ylabel="Population", legend=:right)


# T - Infected population
p_t1 = plot(T1_data.Date[1:horizon], [solsim[4,:], ymeas_u40.T_u40[1:horizon]],
     label=["Simulated Curve" "Measured Data"], title="T - u40", ylabel="Population", legend=:topright)
p_t2 = plot(T1_data.Date[1:horizon], [solsim[10,:], ymeas_mid.T_mid[1:horizon]],
     title="T - 40-60", ylabel="Population", legend =false)
p_t3 = plot(T1_data.Date[1:horizon], [solsim[16,:], ymeas_old.T_old[1:horizon]],
     title="T - 60-80", ylabel="Population", legend =false)
p_t4 = plot(T1_data.Date[1:horizon], [solsim[22,:], ymeas_ger.T_ger[1:horizon]],
     title="T - 80+", ylabel="Population", legend =false)
t_plot = plot(p_t1, p_t2, p_t3, p_t4,size=(900, 600))

T_sim_plot = plot(T1_data.Date[1:horizon], [solsim[4,:],solsim[10,:],solsim[16,:],solsim[22,:]],  label=["Under 40" "Middle-aged (40-60)" "Senior (60-80)" "Geriatric (80+)"], title="T - Threatened", ylabel="Population", legend=:right)
hline!(T_sim_plot, [7e-4], color=:black, linestyle=:dash, label="ICUs Threshold")

# H - Infected population
p_h1 = plot(T1_data.Date[1:horizon], [solsim[5,:], ymeas_u40.H_u40[1:horizon]],
     label=["Simulated Curve" "Measured Data"], title="H - u40", ylabel="Population", legend=:topright)
p_h2 = plot(T1_data.Date[1:horizon], [solsim[11,:], ymeas_mid.H_mid[1:horizon]],
     title="H - 40-60", ylabel="Population", legend =false)
p_h3 = plot(T1_data.Date[1:horizon], [solsim[17,:], ymeas_old.H_old[1:horizon]],
     title="H - 60-80", ylabel="Population", legend =false)
p_h4 = plot(T1_data.Date[1:horizon], [solsim[23,:], ymeas_ger.H_ger[1:horizon]],
     title="H - 80+", ylabel="Population", legend =false)

h_plot = plot(p_h1, p_h2, p_h3, p_h4,size=(900, 600))

# E - Infected population
p_e1 = plot(T1_data.Date[1:horizon], [solsim[6,:], ymeas_u40.E_u40[1:horizon]],
     label=["Simulated Curve" "Measured Data"], title="E - u40", ylabel="Population", legend=:topright)
p_e2 = plot(T1_data.Date[1:horizon], [solsim[12,:], ymeas_mid.E_mid[1:horizon]],
     title="E - 40-60", ylabel="Population", legend =false)
p_e3 = plot(T1_data.Date[1:horizon], [solsim[18,:], ymeas_old.E_old[1:horizon]],
     title="E - 60-80", ylabel="Population", legend =false)
p_e4 = plot(T1_data.Date[1:horizon], [solsim[24,:], ymeas_ger.E_ger[1:horizon]],
     title="E - 80+", ylabel="Population", legend =false)

e_plot = plot(p_e1, p_e2, p_e3, p_e4,size=(900, 600))