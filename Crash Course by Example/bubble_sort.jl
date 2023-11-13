# Bubble Sort 
# This algorithm takes an input array, indexed 1, ..., n, then sort 
# the elements smallest to largest by allowing the larger elements, or "Bubbles",
# to percolate up. 

function bubbleSort!(a)
    n = length(a)
    for i in 1:n-1
        for j in 1:n-i
            if a[j] > a[j+1]
                a[j], a[j+1] = a[j+1], a[j]
            end
        end 
    end 
    return a 
end 

data = [65, 51, 32, 12, 23, 84, 68, 1]

bubbleSort!(data) 
