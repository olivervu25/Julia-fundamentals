module MFABC

export multi_fidelity_abc, Particle

# Define a struct for storing the parameters and weights of a particle
struct Particle
    params::Dict
    weight::Float64
end

# Function to calculate Euclidean distance between two vectors
function euclidean_distance(v1, v2)
    return sqrt(sum((v1 - v2).^2))
end

# Multi-fidelity ABC algorithm
function multi_fidelity_abc(observed_data_lofi, observed_data_hifi, prior_sampler, simulate_lofi, simulate_hifi, epsilon_lofi, epsilon_hifi, eta1, eta2, N)
    particles = Vector{Particle}(undef, N)

    for i = 1:N
        theta = prior_sampler()
        U = rand()

        # Simulate from the low-fidelity model
        simulated_lofi = simulate_lofi(theta)
        weight_tilde = euclidean_distance(simulated_lofi, observed_data_lofi) < epsilon_lofi ? 1.0 : 0.0
        weight = weight_tilde  # Initial weight set to the low-fidelity weight

        # Update the continuation probability
        eta = eta1 * weight_tilde + eta2 * (1 - weight_tilde)

        # Simulate from the high-fidelity model if U < eta
        if U < eta
            simulated_hifi = simulate_hifi(theta)
            weight = euclidean_distance(simulated_hifi, observed_data_hifi) < epsilon_hifi ? 1.0 : 0.0
            # Update weight using the high-fidelity outcome
            weight = weight_tilde + (weight - weight_tilde) / eta
        end

        # Store the particle
        particles[i] = Particle(theta, weight)
    end

    # Calculate weighted average of the parameter function F(theta)
    weighted_sum_F = sum([particle.weight * F(particle.params) for particle in particles])
    weights_sum = sum([particle.weight for particle in particles])
    mu_ABC = weighted_sum_F / weights_sum

    return particles, mu_ABC
end

end
