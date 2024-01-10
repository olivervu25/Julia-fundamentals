using Statistics

function calculate_parameter_means(data)
    num_params = length(data[1]) # Assuming all sub-arrays have the same length
    means = zeros(num_params)

    for i in 1:num_params
        param_values = [data[j][i] for j in 1:length(data)]
        means[i] = mean(param_values)
    end

    return means
end

data = [
    [13.0, 28.0, 53.0, 23.0, 6.0, 31.0],
    [55.0, 62.0, 142.0, 23.0, 23.0, 30.0],
    [75.0, 112.0, 174.0, 42.0, 27.0, 92.0],
    [68.0, 119.0, 205.0, 47.0, 72.0, 135.0],
    [63.0, 125.0, 188.0, 64.0, 106.0, 151.0]
]

parameter_means = calculate_parameter_means(data)

print(parameter_means)