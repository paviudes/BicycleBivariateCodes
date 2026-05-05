function rank_mod2(mat::Matrix{Int})
    """
    Compute the rank of a matrix over GF(2) using Nemo's matrix space functionality.
    """
    field_f2 = Nemo.GF(2)
    nrows, ncols = size(mat)
    ms = matrix_space(field_f2, nrows, ncols)
    mat_rank = rank(ms(mat))
    return mat_rank
end

function nullspace_mod2(mat::Matrix{Int})
    """
    Compute the nullspace of the vector space spanned by the rows of A over GF(2).
    Returns a matrix whose rows form a basis of the nullspace.
    """
    field_f2 = Nemo.GF(2)
    nrows, ncols = size(mat)

    # 2. Create the matrix space and a specific matrix
    ms = matrix_space(field_f2, nrows, ncols)
    mat_field = ms(mat)
    # Compute the nullspace of A_f using Nemo's nullspace function.
    _, N = nullspace(mat_field)

    # Convert to Julia Int matrix
    M = Matrix(N)
    M_int = map(x -> isone(x) ? 1 : 0, M)

    return transpose(M_int)  # rows = basis vectors
end

function compute_logicals(stabilizers_primal::Matrix{Int}, stabilizers_dual::Matrix{Int})::Matrix{Int}
    """
    Given the parity-check matrix (HX or HZ) of a CSS code, compute the logical operators (LZ or LX) of the code.
    """
    # Compute the nullspace of HX to get all operators that commute with the X-stabilizers, i.e. NZ = span(HZ | LZ)
    normalizers = nullspace_mod2(stabilizers_primal)
    nb_normalizers = size(normalizers, 1)

    # println("$(nb_normalizers) operators in N(S)\n", "Normalizers:\n", normalizers, "\n")

    # Greedily extend span(HZ) with normalizer basis vectors.
    # Each accepted candidate is folded into the running span so later candidates are tested
    # against span(HZ ∪ accepted_logicals), preventing multiple representatives of the same
    # logical coset from being counted separately.
    logical_operators = []
    current_span = copy(stabilizers_dual)
    for i in 1:nb_normalizers
        candidate = normalizers[i, :]
        augmented_matrix = vcat(current_span, candidate')
        if rank_mod2(augmented_matrix) > rank_mod2(current_span)
            push!(logical_operators, candidate)
            current_span = copy(augmented_matrix)
        end
    end
    # println("Identified $(length(logical_operators)) logical operators.\n", "Logical operators:\n", logical_operators)
    logical_operators = hcat(logical_operators...)'
    return logical_operators
end

function canonicalize_logicals(logical_X::Matrix{Int}, logical_Z::Matrix{Int})
    """
    Given logical X and Z operators, transform them into a canonical form where the i-th logical X operator anti-commutes only with the i-th logical Z operator and commutes with all others.
     - Compute the symplectic Gram matrix P[i,j] = LX_i · LZ_j mod 2.
     - Find a transformation T such that P @ T^T = I, which can be done by taking T = (P^T)^{-1} over GF(2).
     - Transform the logical Z operators by new_LZ = T * logical_Z mod 2 to achieve the canonical pairing.
     - Return the original logical X operators and the transformed logical Z operators.
    """
    # Symplectic Gram matrix: P[i,j] = LX_i · LZ_j mod 2.
    # The independent X- and Z-logical bases are chosen separately, so P is generally not I.
    # Transform LZ by T = (P^T)^{-1} over GF(2) so that LX @ new_LZ^T = P @ T^T = P @ P^{-1} = I,
    # giving the canonical pairing LX_i anticommutes with LZ_i and commutes with LZ_{j≠i}.
    k = size(logical_X, 1)
    P = mod.(logical_X * logical_Z', 2)
    field_f2 = Nemo.GF(2)
    ms = matrix_space(field_f2, k, k)
    T = map(x -> isone(x) ? 1 : 0, Matrix(inv(transpose(ms(P)))))
    new_LZ = mod.(T * logical_Z, 2)
    return (logical_X, new_LZ)
end

function compute_logical_operators(stabilizers_X::Matrix{Int}, stabilizers_Z::Matrix{Int})
    """
    Given the parity-check matrices HX and HZ of a CSS code, compute the logical operators of the code.
    The logical operators can be found by computing the nullspace of the combined parity-check matrix [HX; HZ] and then identifying which of those operators commute with all stabilizers but are not themselves in the stabilizer group.
    
    Given the Null spaces:
    NX := nullspace_mod2(HX)
    NZ := nullspace_mod2(HZ)

    we know that
    NX = span(HZ | LZ) where LZ are the logical Z operators,
    NZ = span(HX | LX) where LX are the logical X operators.

    To isolate the logical operators LZ from NX, we should find the vectors in NX that are not in the span of HZ.
    This can be done by adding each row of NX to HZ and checking if its rank increases. If it does, then that row corresponds to a logical operator.

    Similarly, to isolate the logical operators LX from NZ, we should find the vectors in NZ that are not in the span of HX.
    """

    logical_Z_operators = compute_logicals(stabilizers_X, stabilizers_Z)
    logical_X_operators = compute_logicals(stabilizers_Z, stabilizers_X)
    canonical_logical_X_operators, canonical_logical_Z_operators = canonicalize_logicals(logical_X_operators, logical_Z_operators)

    return (canonical_logical_X_operators, canonical_logical_Z_operators)
end