# Monte Carlo Simulation 

# Esimate pi 

using Random, LinearAlgebra, Plots; pyplot() 

Random.seed!() 

N = 10^5

data = [[rand(), rand()] for _ in 1:N]

indata = filter((x) -> (norm(x) <= 1), data)
outdata = filter((x) -> (nor(x) > 1), data)

piApprox = 4*length(indata)/N 

println("Pi Estimate: ", piApprox) 


