# Load of the SIDTTHE data from matlab
using MAT
using DataFrames
using Dates
using OrdinaryDiffEq
using Optimization, OptimizationPolyalgorithms, OptimizationOptimJL
using SciMLSensitivity, Zygote, ForwardDiff
using Plots, LaTeXStrings
using Statistics, StatsBase 
using OptimizationMOI, Ipopt

Npop = 59240329  # Total Population of Italy
N = 399          # MPC horizon length
N_mhe = 21       # Estimation horizon (3 weeks)
Ts = 1

# Create the date vector 
start_date = Date(2020, 8, 31) # January 1, 2023
end_date = Date(2021, 10, 3)  # January 10, 2023
date = collect(start_date:end_date)

SIDTTHEfile = matopen("SIDTTHE_data_DEF.mat")
SIDTTHE_data = read(SIDTTHEfile, "SIDTTHE_data")
close(SIDTTHEfile)

I_data = SIDTTHE_data[1, 1]["data"] / Npop
D_data = SIDTTHE_data[2, 1]["data"] / Npop
T_data = (SIDTTHE_data[3, 1]["data"] + SIDTTHE_data[4, 1]["data"]) / Npop
H_data= SIDTTHE_data[5, 1]["data"] / Npop
E_data = SIDTTHE_data[6, 1]["data"] / Npop
S_data = ones(1,399) - (I_data + D_data + T_data + H_data + E_data)

# Gather all the data in a single data matrix
measure_mat = vcat(S_data, I_data, D_data, T_data, H_data, E_data)

function SIDTHEfun!(du, u, p, t)
    s, i, d, t, h, e = u 
    α, δ, γ, σ, τ = p
    λ = 0.1     # coefficient Lambda (fixed)
    du[1] = -s * i * α                         # dS/dt
    du[2] = s * i * α - (γ + λ) * i            # dI/dt
    du[3] = i * γ - d * (λ + δ)                # dD/dt
    du[4] = δ * d - (σ + τ) * t                # dT/dt
    du[5] = (i + d) * λ + t * σ                # dH/dt
    du[6] = t * τ                              # dE/dt
end

u0 = [S_data[1]; I_data[1]; D_data[1]; T_data[1]; H_data[1]; E_data[1]]
α0 = 0.25 
γ0 = 0.12
δ0 = 0.01
σ0 = 0.02
τ0 = 0.02
p0 = [α0, γ0, δ0, σ0, τ0]

# First simulation with initial condition parameters
tspan = (1.0, N_mhe)
tsteps = collect(range(1, stop = N_mhe, length = N_mhe))
prob = ODEProblem(SIDTHEfun!, u0, tspan, p0)

function simulate(p)
    newprob = remake(prob, u0 = u0, p = p)
    solsim = solve(newprob, Tsit5(), saveat = tsteps)
    return reduce(hcat,solsim.u)
end

function loss(p)
    sol = simulate(p)
    # loss = sum(abs2, ( (sol .- ymeas) ./ maximum.(eachrow(ymeas)) )) #+ sum(abs2, p - p_tilde)
    loss = sum(abs2, ( (sol .- ymeas) ./ std.(eachrow(measure_mat)) ))
end

p_lb = [0.05, 0.005, 1e-4, 1e-4, 1e-4] # lower bounds on parameters
p_ub = [0.8, 0.6, 0.6, 0.5, 0.5] # upper bound on parameters

states_dyn = Array{Float64}(undef, 6, 1)
params_dyn = Array{Float64}(undef, 5, 1)

for k in N_mhe:N-N_mhe
#for k in N_mhe:1:N_mhe+10 # Test for three consecutive time windows
    global u0
    global ymeas = measure_mat[:, k-N_mhe+1:k] #available measurement for each time window 

    if k == N_mhe
        global p_tilde = p0
    else 
        u0 = reshape(opti_odesol[:,2],6,1)
        global p_tilde = optires.u
    end

    adtype = Optimization.AutoForwardDiff()
    optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
    optprob = Optimization.OptimizationProblem(optf, p0, lb=p_lb, ub=p_ub) # Initial guess for parameter estimation
    optires = Optimization.solve(optprob, Ipopt.Optimizer(), maxiters = 5000) # Ipopot Optimizer
    opti_odesol = solve(remake(prob, u0=u0, p = optires.u), Tsit5(), saveat = tsteps)
    opti_odesol = reduce(hcat,opti_odesol.u)

    # Save of the optimization values
    params_dyn = hcat( params_dyn, reshape(optires.u, 5, 1) )
    states_dyn = hcat( states_dyn, reshape(opti_odesol[:,end],6,1) )

    # update on the initial conditions
    p0 = reshape(optires.u, 5, 1)
    u0 = reshape(opti_odesol[:,2],6,1)

end

# --- Results Plot ---
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=true) 
scalefontsizes(1)

plotS = plot(date[N_mhe:N-N_mhe], measure_mat[1,N_mhe:N-N_mhe],label="Real Data")
plotS = plot!(date[N_mhe:N-N_mhe], states_dyn[1,2:end],label="Sim Opti Params")
plotS = plot!(ylabel="Population", title="S - Susceptible population")

plotI = plot(date[N_mhe:N-N_mhe], measure_mat[2,N_mhe:N-N_mhe],label="Real Data")
plotI = plot!(date[N_mhe:N-N_mhe], states_dyn[2,2:end],label="Sim Opti Params")
plotI = plot!(ylabel="Population", title="I - Infected population")

plotD = plot(date[N_mhe:N-N_mhe], measure_mat[3,N_mhe:N-N_mhe],label="Real Data")
plotD = plot!(date[N_mhe:N-N_mhe], states_dyn[3,2:end],label="Sim Opti Params")
plotD = plot!(ylabel="Population", title="D - Detected population")

plotT = plot(date[N_mhe:N-N_mhe], measure_mat[4,N_mhe:N-N_mhe],label="Real Data")
plotT = plot!(date[N_mhe:N-N_mhe], states_dyn[4,2:end],label="Sim Opti Params")
plotT = plot!(ylabel="Population", title="T - Threaten population")

plot_α = plot(date[N_mhe:N-N_mhe], params_dyn[1,2:end],label="Real Data")
plot_α = plot!(ylabel="Population", title="α parameter")
