using Pkg 
Pkg.add("DifferentialEquations") 
Pkg.add("Distributions")
using Distributions
using DifferentialEquations

alpha0 = 1
beta = 5 
alpha = 1000
n = 2
Kh = 20
Tfinal = 10

params = (alpha0 = alpha0, beta = beta, alpha = alpha, n = n, Kh = Kh)

#Initial state: m1, m2, m3, p1, p2, p3
initial_state = [0, 0, 0, 40, 20, 60]  

function repressilator(du, u, p, t)
    m1, m2, m3, p1, p2, p3 = u
    alpha0, beta, alpha, n, Kh = p

    # Repression functions
    repression1 = (Kh^n) / (Kh^n + p3^n)
    repression2 = (Kh^n) / (Kh^n + p1^n)
    repression3 = (Kh^n) / (Kh^n + p2^n)

    # Transcription and degradation of m1, m2, m3
    du[1] = alpha0 + alpha * repression1 - m1
    du[2] = alpha0 + alpha * repression2 - m2
    du[3] = alpha0 + alpha * repression3 - m3

    # Translation and degradation of p1, p2, p3
    du[4] = beta * m1 - beta * p1
    du[5] = beta * m2 - beta * p2
    du[6] = beta * m3 - beta * p3
end


function propensities_and_stoichiometric(state)
    m1, m2, m3, p1, p2, p3 = state

    a = [
        alpha0 + alpha * (Kh^n) / (Kh^n + p3^n),  # Transcription of m1
        m1,                                      # Degradation of m1
        alpha0 + alpha * (Kh^n) / (Kh^n + p1^n),  # Transcription of m2
        m2,                                      # Degradation of m2
        alpha0 + alpha * (Kh^n) / (Kh^n + p2^n),  # Transcription of m3
        m3,                                      # Degradation of m3
        beta * m1,                               # Translation of p1
        beta * p1,                               # Degradation of p1
        beta * m2,                               # Translation of p2
        beta * p2,                               # Degradation of p2
        beta * m3,                               # Translation of p3
        beta * p3                                # Degradation of p3
    ]

    stoichiometric = [
        [1, 0, 0, 0, 0, 0],   # Transcription of m1
        [-1, 0, 0, 0, 0, 0],  # Degradation of m1
        [0, 1, 0, 0, 0, 0],   # Transcription of m2
        [0, -1, 0, 0, 0, 0],  # Degradation of m2
        [0, 0, 1, 0, 0, 0],   # Transcription of m3
        [0, 0, -1, 0, 0, 0],  # Degradation of m3
        [0, 0, 0, 1, 0, 0],   # Translation of p1
        [0, 0, 0, -1, 0, 0],  # Degradation of p1
        [0, 0, 0, 0, 1, 0],   # Translation of p2
        [0, 0, 0, 0, -1, 0],  # Degradation of p2
        [0, 0, 0, 0, 0, 1],   # Translation of p3
        [0, 0, 0, 0, 0, -1]   # Degradation of p3
    ]

    return a, stoichiometric
end


# Gillespie SSA algorithm
function gillespie_ssa(state, Tfinal)
    t = 0.0
    times = [t]
    states = [copy(state)]

    while t < Tfinal
        a, updates = propensities_and_stoichiometric(state)
        a0 = sum(a)

        if a0 == 0.0
            break
        end

        tau = rand(Exponential(1/a0))
        t += tau

        if t > Tfinal
            break
        end

        r = rand() * a0
        a_cumsum = cumsum(a)
        reaction_index = findfirst(x -> x > r, a_cumsum)

        state += updates[reaction_index]

        push!(times, t)
        push!(states, copy(state))
    end

    return times, states
end

# Run the simulation
times, states = gillespie_ssa(initial_state, Tfinal)
using Plots

# Assuming `times` and `states` are the outputs from your `gillespie_ssa` function
function plot_gillespie(times, states)
    # Convert the states array to a matrix for easier plotting
    state_matrix = hcat(states...)

    # Create a plot with a title and labels
    p = plot(times, state_matrix', 
             label=["m1" "m2" "m3" "p1" "p2" "p3"], 
             xlabel="Time", 
             ylabel="Concentration", 
             title="Repressilator Model Simulation",
             lw=2)

    return p
end

# Call the plotting function with your simulation results
p = plot_gillespie(times, states)
display(p)
