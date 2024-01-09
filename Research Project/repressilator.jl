module RepressilatorABC

include("mfabc.jl") # Include the definitions from mfabc.jl
include("simulation_method.jl")

using .MFABC         # Use the MFABC module
using .StochasticSimulationModels

using LinearAlgebra
using Distributions

# Define the stoichiometric matrix and initial conditions
const nu = vcat(hcat(I,-I,zeros(3,3),zeros(3,3)), hcat(zeros(3,3),zeros(3,3),I,-I))
const x0 = [0.0, 0.0, 0.0, 40.0, 20.0, 60.0]
const T = 10.0  # Simulation time
const k_nominal = (1.0, 2.0, 5.0, 1000.0, 20.0)


# Define the prior function
function prior()
    k2 = rand(Uniform(1, 4))
    k5 = rand(Uniform(10, 30))
    k_nominal = (1.0, 2.0, 5.0, 1000.0, 20.0)
    return (k_nominal[1], k2, k_nominal[3], k_nominal[4], k5)
end

function repressilator_propensity(state, params)
    m1, m2, m3, p1, p2, p3 = state
    alpha0, beta, alpha, n, Kh = params

    v = zeros(12)
    repression = (Kh^n) ./ (Kh^n .+ [p3^n, p1^n, p2^n])
    v[1:3] = alpha0 .+ alpha .* repression
    v[4:6] = [m1, m2, m3]
    v[7:9] = beta .* [m1, m2, m3]
    v[10:12] = beta .* [p1, p2, p3]

    return v
end

# Instantiate the GillespieModel and TauLeapModel for the repressilator
gillespie_model = GillespieModel(nu, repressilator_propensity, x0, k_nominal, T)

# tau_leap_model = TauLeapModel(nu, repressilator_propensity, x0, k_nominal, T, 0.1)
tau_leap_model = GillespieModel(nu, repressilator_propensity, x0, k_nominal, T)


# Define the low-fidelity and high-fidelity simulation functions
function simulate_lofi(params)
    #simulated_data = simulate_tau_leap(tau_leap_model)
    simulated_data = gillespie_ssa(tau_leap_model)
    return repressilator_summary_statistics(simulated_data)
end

function simulate_hifi(params)
    simulated_data = gillespie_ssa(gillespie_model)
    return repressilator_summary_statistics(simulated_data)
end

# Define the observed data for low-fidelity and high-fidelity models
observed_data_lofi = simulate_lofi(prior())  
observed_data_hifi = simulate_hifi(prior())

# Define epsilon thresholds for low-fidelity and high-fidelity models
epsilon_lofi = 50  # Placeholder value
epsilon_hifi = 50  # Placeholder value

# Define continuation probabilities eta1 and eta2
eta1 = 0.25  # Placeholder value
eta2 = 0.12  # Placeholder value

# Apply the MFABC algorithm
N = 100  # Number of particles to simulate
mfabc_results = multi_fidelity_abc(
    observed_data_lofi, observed_data_hifi, prior, simulate_lofi, simulate_hifi, 
    epsilon_lofi, epsilon_hifi, eta1, eta2, N
)


end 