function verify_canonical_commutations(stabilizers_X::Matrix{Int}, stabilizers_Z::Matrix{Int}, logical_X_operators::Matrix{Int}, logical_Z_operators::Matrix{Int})
    """
    Test if the logical X and Z operators computed from the stabilizer matrices commute correctly.
    - all X stabilizers and logical Z operators should commute (i.e. their symplectic product should be 0 mod 2)
    - all Z stabilizers and logical X operators should commute (i.e. their symplectic product should be 0 mod 2)
    - each logical X operator should anti-commute with its corresponding logical Z operator (i.e. their symplectic product should be 1 mod 2)
    """
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