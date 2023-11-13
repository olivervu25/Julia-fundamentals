# Monte Carlo Simulation 

# Esimate pi 

using Random, LinearAlgebra, Plots; pyplot() 

Random.seed!() 

N = 10^5

data = [[rand(), rand()] for _ in 1:N]

indata = filter((x) -> (norm(x) <= 1), data)
outdata = filter((x) -> (norm(x) > 1), data)

piApprox = 4*length(indata)/N 

println("Pi Estimate: ", piApprox) 

# first.(indata), last.(indata) are used to extract the x and y coordinates of points within the indata
# c:=blue sets the color of these points to blue 
# ms=1 and msw=0 set the marker size and marker stroke width, respectively, making points smaller and without an outline 
scatter(first.(indata),last.(indata), c=:blue, ms=1, msw=0)
scatter!(first.(outdata),last.(outdata), c=:red, ms=1, msw=0, xlims=(0,1), ylims=(0,1), legend=:none, ratio=:equal)