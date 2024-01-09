using DifferentialEquations
using Distributions  # For the Poisson distribution
using Random

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

#Using Poisson Tau-Leap with tau=1 for Low-fidelity 
function simulate_low_fidelity(initial_state, Tfinal, params)
    t = 0.0
    times = [t]
    states = [copy(initial_state)]

    while t < Tfinal
        a, updates = propensities_and_stoichiometric(states[end], params)
        a0 = sum(a)

        if a0 == 0.0
            break
        end

        # Determine a suitable time step `tau` (you may need to calculate this)
        tau = 0.5

        # Estimate the number of events for each reaction using the Poisson distribution
        num_events = [rand(Poisson(a[i] * tau)) for i in 1:length(a)]

        # Update the state according to the number of events
        for (reaction_idx, events) in enumerate(num_events)
            state_changes = updates[reaction_idx] .* events
            states[end] .+= state_changes
        end

        # Update time and record the new state
        t += tau
        push!(times, t)
        push!(states, copy(states[end]))
    end

    return times, states
end

#Using SSA Implementation for high fidelity 
function simulate_high_fidelity(initial_state, Tfinal, params)
    t = 0.0
    times = [t]
    states = [copy(initial_state)]

    while t < Tfinal
        a, updates = propensities_and_stoichiometric(states[end], params)
        a0 = sum(a)

        if a0 == 0.0
            break
        end

        # Gillespie SSA chooses the time to next reaction using an exponential distribution
        tau = rand(Exponential(1/a0))

        # It also chooses which reaction will occur based on the propensities
        r = rand() * a0
        a_cumsum = cumsum(a)
        reaction_index = findfirst(x -> x > r, a_cumsum)

        # Update the state based on the chosen reaction
        state_changes = updates[reaction_index]
        states[end] .+= state_changes

        # Update time and record the new state
        t += tau
        if t > Tfinal
            break
        end
        push!(times, t)
        push!(states, copy(states[end]))
    end

    return times, states
end

#Prior 
function prior()
    return Dict(
        :alpha0 => 1, 
        :beta => 5,
        :alpha => 1000,
        :n => rand(Uniform(1, 4)),  # n ~ U(1, 4)
        :Kh => rand(Uniform(10, 30)) # Kh ~ U(10, 30)
    )
end


# Multi-fidelity ABC function (full version)
function mf_abc(observed_data_lofi, observed_data_hifi, Tfinal, N, tau_lofi, tau_hifi, epsilon_lofi, epsilon_hifi, eta1, eta2)
    weights = zeros(N)
    theta_list = [Dict() for _ in 1:N]
    distance_lofi = zeros(N)
    distance_hifi = zeros(N)

    # Tracking acceptance
    accepted_lofi = Bool[]
    accepted_hifi = Bool[]
    false_positive = Bool[]
    false_negative = Bool[]

    for i in 1:N
        theta_i = prior()
        U = rand()

        # Low-fidelity model simulation
        simulated_data_lofi = simulate_low_fidelity(state, Tfinal, theta_i)
        distance_lofi[i] = distance(simulated_data_lofi, observed_data_lofi)
        w_tilde = distance_lofi[i] < epsilon_lofi
        w_i = w_tilde ? 1.0 : 0.0
        push!(accepted_lofi, w_tilde)

        # Continuation with high-fidelity model
        eta = eta1 * w_i + eta2 * (1 - w_i)
        if U < eta
            simulated_data_hifi = simulate_high_fidelity(state, Tfinal, theta_i, tau_hifi)
            distance_hifi[i] = distance(simulated_data_hifi, observed_data_hifi)
            w = distance_hifi[i] < epsilon_hifi
            w_i *= w ? 1.0 : (1 - eta1) / eta2
            push!(accepted_hifi, w)
            push!(false_positive, w_tilde && !w)
            push!(false_negative, false)  # Cannot have false negatives if high-fidelity is evaluated
        else
            distance_hifi[i] = NaN  # High-fidelity was not evaluated
            push!(accepted_hifi, false)
            push!(false_positive, false)
            push!(false_negative, !w_tilde)  # Potential false negative if rejected by low-fidelity
        end

        weights[i] = w_i
        theta_list[i] = theta_i
    end

    # Compute weighted parameters
    weighted_params = sum(weights .* map(theta -> theta[:n], theta_list)) / sum(weights)

    return weights, theta_list, weighted_params, distance_lofi, distance_hifi, accepted_lofi, accepted_hifi, false_positive, false_negative
end


N = 100  # Number of samples in ABC
tau_lofi = 0.5  # Time step for low-fidelity simulation
tau_hifi = 0.1  # Time step for high-fidelity simulation
epsilon_lofi = 200  # Epsilon threshold for low-fidelity model
epsilon_hifi = 200  # Epsilon threshold for high-fidelity model
eta1 = 0.25  # Continuation probability for lo-fi to hi-fi
eta2 = 0.12  # Continuation probability for hi-fi model check

results = mf_abc(observed_data_lofi, observed_data_hifi, Tfinal, N, tau_lofi, tau_hifi, epsilon





