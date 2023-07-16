#==Balancing chemical equations programmatically================================

To balance a chemical equation using matrix operations, we first need to represent the equation as a system of linear equations. 
We can do this by creating a matrix that represents the coefficients of the elements in each compound. 
For example, consider the equation `CH4 + O2 -> CO2 + H2O`. We can represent this equation as a system of linear equations:

```
C: 1x + 0y - 1z = 0
H: 4x + 0y - 2w = 0
O: 2y - 2z - 1w = 0
```

where `x`, `y`, `z`, and `w` are the coefficients of `CH4`, `O2`, `CO2`, and `H2O`, respectively. 
The first equation represents the conservation of carbon atoms, 
the second equation represents the conservation of hydrogen atoms, 
and the third equation represents the conservation of oxygen atoms.

We can represent this system of linear equations as a matrix:

```
1 0 -1 0
4 0 0 -2
0 2 -2 -1
```

Each row corresponds to an equation and each column corresponds to a compound. 
The first column corresponds to `CH4`, the second column corresponds to `O2`, 
the third column corresponds to `CO2`, and the fourth column corresponds to `H2O`.

To solve for the coefficients of the balanced equation, we need to find the nullspace of the matrix. 
The nullspace of a matrix is the set of all vectors that, when multiplied by the matrix, result in a zero vector. 
In other words, we're looking for the values of `x`, `y`, `z`, and `w` that satisfy the system of linear equations.

===============================================================================#

using LinearAlgebra

function balance(equation::String)
    # Validate equation
    if !('=' in equation)
        return "Invalid equation"
    end

    invalidChars::Array{SubString{String}} = split("`~!@#\$%^&*_-[]\\{}|;':\",./<>?")

    for sign in invalidChars
        if occursin(sign, equation)
            return "Invalid sign: " * sign
        end
    end

    # Add elements to Set
    elements::Set = Set()

    for i in eachindex(equation)
        if i < length(equation)
            if isuppercase(equation[i])
                if isuppercase(equation[i + 1]) || tryparse(Int8, string(equation[i + 1])) !== nothing
                    push!(elements, string(equation[i]))
                elseif islowercase(equation[i + 1])
                    push!(elements, equation[i] * equation[i + 1])
                end
            end
        else
            if isuppercase(equation[i])
                push!(elements, string(equation[i]))
            end
        end
    end
    
    # Count elements
    equalSign::Int = findfirst("=", equation)[1]
    substances::Array = split(equation, r"[\+\=]")
    linearEqs::Matrix = zeros(Int64, length(elements), length(substances))
    i = 1

    for element in elements
        for x in eachindex(substances)
            occurences = findall(element, substances[x])
            count = 0

            for index in occurences
                leftParens = 0
                skip = 0
                subCount = 1

                for j in (index[end] + 1): length(substances[x])
                    if tryparse(Int128, string(substances[x][j])) === nothing && j > index[end] + 1
                        subCount = tryparse(Int128, substances[x][(index[end] + 1):(j - 1)])
                        break
                    elseif tryparse(Int128, string(substances[x][j])) === nothing
                        break
                    end
                end

                for j in index[1]:-1:1
                    if substances[x][j] == '('
                        leftParens += 1
                    end
                end

                for j in index[end]:length(substances[x])
                    if substances[x][j] == ')'
                        if skip == 0
                            for k in (j + 1): length(substances[x])
                                if tryparse(Int128, string(substances[x][k])) === nothing && k > j + 1
                                    subCount = subCount * tryparse(Int128, substances[x][(j + 1):(k - 1)])
                                end
                            end
                        else
                            skip -= 1
                        end
                    elseif substances[x][j] == '('
                        skip += 1
                    end
                end

                count += subCount
            end

            linearEqs[i, x] = count * (equalSign > findall(substances[x], equation)[1][1] ? 1 : -1)
        end

        i += 1
    end

    # Nullspace and round numbers
    coeffs = nullspace(linearEqs)
    minCoeff = min(coeffs...)
    for i in eachindex(coeffs)
        coeffs[i] = round(coeffs[i] * 1 / minCoeff, sigdigits=3)
    end

    # Return balanced equation
    for i in eachindex(substances)
        substances[i] = string(coeffs[i]) * strip(substances[i])
    end

    return join(substances, " : ")
end

println("** [Chemical Equations Balancer] **")
println("###################################")
println(" - Spaces required around \"+\"")
println(" - Two-letter elements only")
println(" - No ions")
println(" - Only round brackets are allowed")
println(" - CH3- => Me; C6H5- => Ph")
println("-----------------------------------")
println("E.g. CH4 + O2 = CO2 + H2O")
println("-----------------------------------")
print("Enter equation: ")
input = readline()
println("___________________________________")
println("+ Result:")
println("___________________________________")
# println(balance(input))

# Examples
# println(balance("CH4 + O2 = CO2 + H2O"))
# println(@time balance("Cu + HNO3 = Cu(NO3)2 + NO + H2O"))
