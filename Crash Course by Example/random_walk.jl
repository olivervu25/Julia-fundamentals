# Random Walk 

using Plots, Random, Measures; pyplot() 

function path(n=5000)
    x, y = 0.0, 0.0
    xDat, yDat = [], [] 
    for _ in 1:n
        flip = rand(1:4)
        if flip == 1
            x += 1
        elseif flip == 2
            y += 1
        elseif flip == 3
            x -= 1
        elseif flip == 4
            y -= 1
        end 
        push!(xDat,x)
        push!(yDat,y)

        println(typeof(flip)) 
    end 
    return xDat, yDat
end 



default(xlabel = "x", ylabel = "y", xlims=(-150,50), ylims=(-250,50))
p1 = plot(path(), c=:blue)
p1 = plot!(path(), c=:red)
p1 = plot!(path(), c=:green)

