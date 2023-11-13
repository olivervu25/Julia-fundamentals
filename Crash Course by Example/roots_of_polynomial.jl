# Roots of a Polynomial

# We can add, update, remove packages in Julia 

# ] add Roots 
# ] status: lists what packages and versions are currently installed 
# ] update: updates existing packages 
# ] remove Roots: removes package from the current Julia build 

using Roots 

function polynomialGenerator(a...)
    n = length(a) - 1
    poly = function(x)
        return sum([a[i+1]*x^i for i in 0:n])
    end 

    return poly 
end 

polynomial = polynomialGenerator(1, 3, -10)

# Function find_zeros() returns the roots of the original polynomial 
zeroVals = find_zeros(polynomial,-10,10)

println("Zeros of the function f(x): ", zeroVals)


