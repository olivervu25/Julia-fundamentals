using DifferentialEquations
using Distributions  # For the Poisson distribution

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

#Using Poisson Tau-Leap for Low-fidelity 
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

#Using Binomial tau-leap for high-fidelity 
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

        tau = 0.5

        # Update the state using the binomial tau-leap method
        for (reaction_idx, update) in enumerate(updates)
            # Find the number of available reactant molecules
            num_reactants = # ... calculate based on the state and update
            # Calculate the probability of reaction
            prob_reaction = 1 - exp(-a[reaction_idx] * tau)
            # Estimate the number of events using the Binomial distribution
            num_events = rand(Binomial(num_reactants, prob_reaction))
            state_changes = update .* num_events
            states[end] .+= state_changes
        end

        # Update time and record the new state
        t += tau
        push!(times, t)
        push!(states, copy(states[end]))
    end

    return times, states
end

