using MAT, CSV
using DataFrames
using Dates
using OrdinaryDiffEq, LinearAlgebra
using Optimization, OptimizationPolyalgorithms, OptimizationOptimJL
using SciMLSensitivity, Zygote, ForwardDiff
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

# Manipulation of raw df 
I_data = vcat(I_raw[1:1, :], I_raw)
T_data = DataFrame( Date = T1_data.Date, u40 = T1_data.u40 + T2_data.u40, mid = T1_data.mid + T2_data.mid,
    old = T1_data.old + T2_data.old, ger = T1_data.ger + T2_data.ger
)

# Contac matrices from POLYMOD study 
c_home =  [  19.7896   7.35764  1.43232   1.3;
             4.73514   5.62234  0.445607  0.44;
             3.26361   2.32613  4.3079    4.00; 
             3.00      2.00     4.00      3.00  ]

c_schl =  [  24.9874    2.216       0.0889567   0.0389567;
             8.94216    1.1419      0.0793523   0.0293523;
             0.666502   0.0661153   0.054174    0.014174;
             0.0566502  0.0361153   0.014174    0.004174 ]

c_work =  [  13.2266     8.39663   0.0329148    0.0229148;
             9.40362     8.72244   0.0357202    0.0157202;
             0.0930825   0.154919  0.000840544  0.00040544;
             0.030825    0.015491  0.00040544   0.00000001 ]

c_othr =  [  44.1768     12.4868   4.61979   4.61979
             11.3875     8.96874   4.04791   4.04791
             6.62855     9.04265   8.95355   8.00000 
             3.62855     6.04265   9.00000   9.95355 ]

C = c_home + c_schl + c_work + c_othr

# Guessed values of lambda, assuming that older people heal slower
λ1 = 1/7; λ2 = 1/12; λ3= 1/15; λ4 = 1/15; 

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
    du[1]  = -s1 * α1 * (C[1, 1] * i1 + C[1, 2] * i2 + C[1, 3] * i3 + C[1, 4] * i4)                   # dS1/dt
    du[2]  = -s1 * α1 * (C[1, 1] * i1 + C[1, 2] * i2 + C[1, 3] * i3 + C[1, 4] * i4) - (γ1 + λ1) * i1  # dI1/dt
    du[3]  = i1 * γ1 - d1 * (λ1 + δ1)                                                                 # dD1/dt
    du[4]  = δ1 * d1 - (σ1 + τ1) * t1                                                                 # dT1/dt
    du[5]  = (i1 + d1) * λ1 + t1 * σ1                                                                 # dH1/dt
    du[6]  = t1 * τ1                                                                                  # dE1/dt

    # Set of equations for 40-60 (middle aged) (group 2)
    du[7]  = -s2 * α2 * (C[2, 1] * i1 + C[2, 2] * i2 + C[2, 3] * i3 + C[2, 4] * i4)                   # dS2/dt
    du[8]  = -s2 * α2 * (C[2, 1] * i1 + C[2, 2] * i2 + C[2, 3] * i3 + C[2, 4] * i4) - (γ2 + λ2) * i2  # dI2/dt
    du[9]  = i2 * γ2 - d2 * (λ2 + δ2)                                                                 # dD2/dt
    du[10] = δ2 * d2 - (σ2 + τ2) * t2                                                                 # dT2/dt
    du[11] = (i2 + d2) * λ2 + t2 * σ2                                                                 # dH2/dt
    du[12] = t2 * τ2                                                                                  # dE2/dt

    # Set of equations for 60-80 (old aged) (group 3)
    du[13] = -s3 * α3 * (C[3, 1] * i1 + C[3, 2] * i2 + C[3, 3] * i3 + C[3, 4] * i4)                   # dS3/dt
    du[14] = -s3 * α3 * (C[3, 1] * i1 + C[3, 2] * i2 + C[3, 3] * i3 + C[3, 4] * i4) - (γ3 + λ3) * i3  # dI3/dt
    du[15] = i3 * γ3 - d3 * (λ3 + δ3)                                                                 # dD3/dt
    du[16] = δ3 * d3 - (σ3 + τ3) * t3                                                                 # dT3/dt
    du[17] = (i3 + d3) * λ3 + t3 * σ3                                                                 # dH3/dt
    du[18] = t3 * τ3                                                                                  # dE3/dt

    # Set of equations for 80+ (group 4)
    du[19] = -s4 * α4 * (C[4, 1] * i1 + C[4, 2] * i2 + C[4, 3] * i3 + C[4, 4] * i4)                   # dS4/dt
    du[20] = -s4 * α4 * (C[4, 1] * i1 + C[4, 2] * i2 + C[4, 3] * i3 + C[4, 4] * i4) - (γ4 + λ4) * i4  # dI4/dt
    du[21] = i4 * γ4 - d4 * (λ4 + δ4)                                                                 # dD4/dt
    du[22] = δ4 * d4 - (σ4 + τ4) * t4                                                                 # dT4/dt
    du[23] = (i4 + d4) * λ4 + t4 * σ4                                                                 # dH4/dt
    du[24] = t4 * τ4                                                                                  # dE4/dt
end

# Initial conditions for the simulation 
u0 = [S_data.u40[1], I_data.u40[1], D_data.u40[1], T_data.u40[1], H_data.u40[1], E_data.u40[1], 
      S_data.mid[1], I_data.mid[1], D_data.mid[1], T_data.mid[1], H_data.mid[1], E_data.mid[1], 
      S_data.old[1], I_data.old[1], D_data.old[1], T_data.old[1], H_data.old[1], E_data.old[1], 
      S_data.ger[1], I_data.ger[1], D_data.ger[1], T_data.ger[1], H_data.ger[1], E_data.ger[1] ]

α01 = 0.1;    α02 = 0.25;  α03 = 0.45; α04 = 0.55
γ01 = 0.1;    γ02 = 0.15;  γ03 = 0.25; γ04 = 0.25
δ01 = 0.005;  δ02 = 0.05;  δ03 = 0.15; δ04 = 0.25
σ01 = 0.5;    σ02 = 0.1;   σ03 = 0.02; σ04 = 0.005
τ01 = 0.0001; τ02 = 0.007; τ03 = 0.25; τ04 = 0.3

p0 = [α01, α02, α03, α04,
      γ01, γ02, γ03, γ04,
      δ01, δ02, δ03, δ04,
      σ01, σ02, σ03, σ04,
      τ01, τ02, τ03, τ04 ]

horizon = 50
tspan = (1.0, horizon)
tsteps = collect(range(1, stop = horizon, length = horizon))
prob = ODEProblem(SIDTHEage!, u0, tspan, p0)
solsim = solve(prob, Tsit5(), saveat = tsteps)
reduce(hcat,solsim.u)
