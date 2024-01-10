module RepressilatorABC

include("mfabc.jl") # Include the definitions from mfabc.jl
include("simulation_method.jl")

using Dierckx

using .MFABC         # Use the MFABC module
using .StochasticSimulationModels
using Plots

using LinearAlgebra
using Distributions

# Define the stoichiometric matrix and initial conditions
const nu = vcat(hcat(I,-I,zeros(3,3),zeros(3,3)), hcat(zeros(3,3),zeros(3,3),I,-I))
const x0 = [0.0, 0.0, 0.0, 40.0, 20.0, 60.0]
const T = 10.0  # Simulation time
const k_nominal = (1.0, 5, 1000, 2, 20.0)


# Define the prior function
function prior()
    n = rand(Uniform(1, 4))
    Kh = rand(Uniform(10, 30))
    #k_nominal = (1.0, 2.0, 5.0, 1000.0, 20.0)
    return (1.0, 5.0, 1000.0, n, Kh) 
end

function repressilator_propensity(state, params)
    m1, m2, m3, p1, p2, p3 = state
    alpha0, beta, alpha, n, Kh = params

    v = zeros(12)
    repression1 = (Kh^n) / (Kh^n + p3^n)
    repression2 = (Kh^n) / (Kh^n + p1^n)
    repression3 =(Kh^n) / (Kh^n + p2^n)
    
    v[1] = alpha0 + alpha * repression1
    v[2] = alpha0 + alpha * repression2
    v[3] = alpha0 + alpha * repression3

    v[4:6] = [m1, m2, m3]
    v[7:9] = beta .* [m1, m2, m3]
    v[10:12] = beta .* [p1, p2, p3]

    return v
end

# Instantiate the GillespieModel and TauLeapModel for the repressilator
gillespie_model = GillespieModel(nu, repressilator_propensity, x0, k_nominal, T)

# tau_leap_model = TauLeapModel(nu, repressilator_propensity, x0, k_nominal, T, 0.1)
tau_leap_model = TauLeapModel(nu, repressilator_propensity, x0, k_nominal, T)

# Assuming `times` and `states` are the outputs from your `gillespie_ssa` function
function plot_gillespie(times, states)
    # Convert the states array to a matrix for easier plotting
    state_matrix = hcat(states...)

    # Create a plot with a title and labels
    p = plot(state_matrix', 
             label=["m1" "m2" "m3" "p1" "p2" "p3"], 
             xlabel="Time", 
             ylabel="Concentration", 
             title="Repressilator Model Simulation",
             lw=2)

    return p
end

# Call the plotting function with your simulation results
#p = plot_gillespie(gillespie_ssa(gillespie_model)[1], gillespie_ssa(gillespie_model)[2])
#display(p)


#p = plot_gillespie(simulate_tau_leap(tau_leap_model)[1], simulate_tau_leap(tau_leap_model)[2])
#display(p)

# Define the low-fidelity and high-fidelity simulation functions
function simulate_lofi(params)
    simulated_data = simulate_tau_leap(tau_leap_model)
    #simulated_data = simulate_tau_leap(tau_leap_model)
    return simulated_data
end

function simulate_hifi(params)
    simulated_data = gillespie_ssa(gillespie_model)
    return simulated_data
end

# Define the observed data for low-fidelity and high-fidelity models
observed_data_lofi = simulate_lofi(prior())  
observed_data_hifi = simulate_hifi(prior())

# Define epsilon thresholds for low-fidelity and high-fidelity models
epsilon_lofi = 1500  # Placeholder value
epsilon_hifi = 7.5  # Placeholder value

# Define continuation probabilities eta1 and eta2
eta1 = 0.25  # Placeholder value
eta2 = 0.12  # Placeholder value

# Apply the MFABC algorithm
N = 5000  # Number of particles to simulate
mfabc_results = multi_fidelity_abc(
    observed_data_lofi, observed_data_hifi, prior, simulate_lofi, simulate_hifi, 
    epsilon_lofi, epsilon_hifi, eta1, eta2, N
)

#print(mfabc_results)

#print(mfabc_results[1]) 

function plot_weights_distribution(particles)
    weights = [p.weight for p in particles]
    histogram(weights, bins=50, title="Weights Distribution", xlabel="Weights", ylabel="Count")
end

# Example usage:
weight_result = []
for particle in mfabc_results
    push!(weight_result,particle.weight) 
end 

p = histogram(weight_result, title="Weights Distribution", xlabel="Weights", ylabel="Count")

display(p)
# Example: Plot weighted means
#plot_weighted_means(mfabc_results)
print("end") 

end 