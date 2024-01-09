module StochasticSimulationModels

export GillespieModel, TauLeapModel, gillespie_ssa, simulate_tau_leap

using Distributions

"""
GillespieModel

A struct representing a model for simulating biochemical reactions using the Gillespie Stochastic 
Simulation Algorithm (SSA). The Gillespie SSA is an exact simulation method that models the 
stochastic behavior of reactions in well-mixed biochemical systems.
"""

Parameters = Tuple{Vararg{Float64}}

struct GillespieModel
    nu::Matrix{Float64} # Stoichiometric matrix
    propensity::Function # Function to calculate propensities
    x0::Array{Float64,1} # Initial conditions
    params::Parameters # Parameters for the model
    T::Float64 # Simulation time
end
    
    
function gillespie_ssa(model::GillespieModel)
    t = 0.0
    times = [t]
    states = [copy(model.x0)]
        
    while t < model.T
        a = model.propensity(states[end], model.params)
        a0 = sum(a)
    
        if a0 == 0.0
            break
        end
    
        tau = rand(Exponential(1/a0))
        t += tau
        if t > model.T
            break
        end
    
        r = rand() * a0
        a_cumsum = cumsum(a)
        reaction_index = findfirst(>(r), a_cumsum)
    
        state_change = model.nu[:, reaction_index]
        states[end] .+= state_change
    
        push!(times, t)
        push!(states, copy(states[end]))
    end
    
    return times, states
    
end


"""
Tau Leap Model 

A struct representing a model for simulating biochemical reactions using the tau-leap method. 
    
The tau-leap method is an approximate stochastic simulation algorithm that allows for larger 
time steps than the exact Gillespie Stochastic Simulation Algorithm (SSA), improving computational 
efficiency at the cost of some accuracy.
"""
 
struct TauLeapModel
    nu::Matrix{Float64}             # Stoichiometric matrix
    propensity::Function            # Function to calculate propensities
    x0::Array{Float64,1}            # Initial conditions
    params::Parameters              # Parameters for the model
    T::Float64                      # Simulation time
    tau::Float64                    # Fixed time step for tau-leap
end

# struct TauLeapModel
#     nu::Matrix{Float64}             # Stoichiometric matrix
#     propensity::Function            # Function to calculate propensities
#     diff_propensity::Function       # Differentiated propensity function
#     x0::Array{Float64,1}            # Initial conditions
#     params::Parameters              # Parameters for the model
#     T::Float64                      # Simulation time
#     tau::Float64                    # Base non-adaptive tau-leap
#     nc::Int                         # Critical reaction threshold
#     epsilon::Float64                # Error parameter for adaptive stepping
# end


function simulate_tau_leap(model::TauLeapModel)
    t = 0.0
    times = [t]
    states = [copy(model.x0)]

    while t < model.T
        a = model.propensity(states[end], model.params)

        # Check for non-negative, non-NaN propensities
        if any(a .< 0) || any(isnan.(a))
            error("Invalid propensity encountered")
        end

        a0 = sum(a)

        if a0 == 0.0
            break
        end

        num_reactions = [rand(Poisson(a_i * model.tau)) for a_i in a]
        state_change = sum(model.nu .* num_reactions, dims=2)
        states[end] .+= state_change[:]

        t += model.tau
        if t > model.T
            break
        end

        push!(times, t)
        push!(states, copy(states[end]))
    end

    return times, states
end

end 