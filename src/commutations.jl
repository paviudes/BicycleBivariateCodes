function verify_canonical_commutations(stabilizers_X::Matrix{Int}, stabilizers_Z::Matrix{Int})
    """
    Test if the logical X and Z operators computed from the stabilizer matrices commute correctly.
    - all X stabilizers and logical Z operators should commute (i.e. their symplectic product should be 0 mod 2)
    - all Z stabilizers and logical X operators should commute (i.e. their symplectic product should be 0 mod 2)
    - each logical X operator should anti-commute with its corresponding logical Z operator (i.e. their symplectic product should be 1 mod 2)
    """
    # println("Stabilizer X matrix shape: $(size(stabilizers_X))\n", "Stabilizer X matrix:\n", stabilizers_X, "\n")
    # println("Stabilizer Z matrix shape: $(size(stabilizers_Z))\n", "Stabilizer Z matrix:\n", stabilizers_Z, "\n")

    # Compute the logical operators using the function.
    logical_X_operators, logical_Z_operators = compute_logical_operators(stabilizers_X, stabilizers_Z)

    # println("Logical X operators shape: $(size(logical_X_operators))\n", "Logical X operators:\n", logical_X_operators, "\n")
    # println("Logical Z operators shape: $(size(logical_Z_operators))\n", "Logical Z operators:\n", logical_Z_operators, "\n")

    # Check commutation between X stabilizers and logical Z operators.
    n_stabilizers_X = size(stabilizers_X, 1)
    n_logical_Z = size(logical_Z_operators, 1)
    for i in 1:n_stabilizers_X
        for j in 1:n_logical_Z
            symplectic_product = mod(sum(stabilizers_X[i, :] .* logical_Z_operators[j, :]), 2)
            if symplectic_product != 0
                println("ERROR: X stabilizer $(i) does not commute with logical Z operator $(j)")
                return false
            end
        end
    end

    # Check commutation between Z stabilizers and logical X operators.
    n_stabilizers_Z = size(stabilizers_Z, 1)
    n_logical_X = size(logical_X_operators, 1)
    for i in 1:n_stabilizers_Z
        for j in 1:n_logical_X
            symplectic_product = mod(sum(stabilizers_Z[i, :] .* logical_X_operators[j, :]), 2)
            if symplectic_product != 0
                println("ERROR: Z stabilizer $(i) does not commute with logical X operator $(j)")
                return false
            end
        end
    end

    # Check anti-commutation between logical X and Z operators.
    for i in 1:n_logical_X
        for j in 1:n_logical_Z
            symplectic_product = mod(sum(logical_X_operators[i, :] .* logical_Z_operators[j, :]), 2)
            if i == j
                if symplectic_product != 1
                    println("ERROR: Logical X operator $(i) does not anti-commute with its corresponding logical Z operator $(j)")
                    return false
                end
            else
                if symplectic_product != 0
                    println("ERROR: Logical X operator $(i) does not commute with logical Z operator $(j) when it should")
                    return false
                end
            end
        end
    end
    
    println("All commutation tests passed!")
    return true
end