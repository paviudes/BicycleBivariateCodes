# ----------------------------------------------------------------------------
# Bivariate Bicycle (BB) quantum LDPC codes from
#   Bravyi, Cross, Gambetta, Maslov, Rall, Yoder (2024),
#   "High-threshold and low-overhead fault-tolerant quantum memory."
#   https://arxiv.org/abs/2308.07915
#
# A BB code is parameterized by integers (l, m) and two polynomials
#   A(x, y), B(x, y) ∈ F_2[x, y] / (x^l - 1, y^m - 1)
# yielding a CSS code with parameters [[2lm, k, d]].
#
# Stabilizer parity-check matrices:
#   H_X = [A | B]      (lm rows × 2lm cols)
#   H_Z = [B^T | A^T]  (lm rows × 2lm cols)
# X-stabilizers act on qubits indexed 1..lm (left block) and lm+1..2lm (right
# block) according to the 1-entries of each row of H_X; similarly for Z.
# ----------------------------------------------------------------------------

struct BBCode
    HX::Matrix{Int}
    HZ::Matrix{Int}
    LX::Matrix{Int}
    LZ::Matrix{Int}
    n::Int
    k::Int
    d::Int
    l::Int
    m::Int
    A_monomials::Vector{Tuple{Int,Int}}
    B_monomials::Vector{Tuple{Int,Int}}
    function BBCode(HX, HZ, LX, LZ, n, k, d, l, m, A_monomials, B_monomials)
        new(HX, HZ, LX, LZ, n, k, d, l, m, A_monomials, B_monomials)
    end
end

# ============================================================================
# BB code construction.
# ============================================================================

function _cyclic_shift(n::Int)::Matrix{Int}
    """
    Build the n × n cyclic shift matrix S, given by
    S[i, j] = 1 if j ≡ i + 1 (mod n)
            = 0 otherwise.
    """
    cyclic = zeros(Int, n, n)
    @inbounds for i in 1:n
        cyclic[i, mod1(i + 1, n)] = 1
    end
    return cyclic
end

function _polynomial_matrix(l::Int, m::Int, monomials::Vector{Tuple{Int,Int}})::Matrix{Int}
    """
    Build the l m × l m matrix corresponding to a bivariate polynomial A(x, y) = ∑_(a,b) x^a y^b.
    x and y are defined as
    x = S_l ⊗ I_m
    y = I_l ⊗ S_m
    where S_n is the n × n cyclic shift matrix given by `_cyclic_shift(n)`.
    Hence:
    x^a y^b = (S_l^a ⊗ I_m) * (I_l ⊗ S_m^b).
            = S_l^a ⊗ S_m^b
    """
    Sl = _cyclic_shift(l)
    Sm = _cyclic_shift(m)
    Il = Matrix{Int}(I, l, l)
    Im = Matrix{Int}(I, m, m)

    # Compute the matrix corresponding to the polynomial A(x, y) = ∑_(a,b) x^a y^b.
    n_qubits = l * m
    poly_matrix = zeros(Int, n_qubits, n_qubits)
    for (a, b) in monomials
        Sla = a == 0 ? Il : Sl^a # matrix power.
        Smb = b == 0 ? Im : Sm^b # matrix power.
        monomial_term = kron(Sla, Smb) # Kronecker product.
        poly_matrix = poly_matrix .+ monomial_term
    end
    
    # Since we are working over F_2, we take the result modulo 2.
    binary_mat = mod.(poly_matrix, 2)
    
    return binary_mat
end

function gen_bb_code(
    l::Int, m::Int,
    A_monomials::Vector{Tuple{Int,Int}},
    B_monomials::Vector{Tuple{Int,Int}};
    distance::Int = -1 # optional distance parameter, since computing the distance is expensive. If not provided, it will be set to -1.
)
    """
    Construct a Bicycle Bivariate (BB) code given by parameters (l, m) and two bivariate polynomials A and B, following the prescription in Section 4 of https://arxiv.org/abs/2308.07915.

    The code is constructed as follows:
    1. Construct a cycling shift matrix S_l of size l × l and S_m of size m × m.
    2. Construct the matrices: x = S_l ⊗ I_m and y = I_l ⊗ S_m.
    3. Construct the matrices A and B corresponding to the input polynomials, as A = ∑_(a,b) x^a y^b and B = ∑_(a,b) x^a y^b.
    4. Construct the parity-check matrices H_X = [A | B] and H_Z = [B^T | A^T].

    The resulting [[n,k,d]] BB code has
    - n = 2 l m physical qubits
    - k = 2 l m - rank(H_X) - rank(H_Z) logical qubits
    - d inferred from Table 3 of https://arxiv.org/abs/2308.07915.

    We will return the parity-check matrices H_X and H_Z, and the [[n, k, d]] parameters of the code
    """
    A = _polynomial_matrix(l, m, A_monomials)
    B = _polynomial_matrix(l, m, B_monomials)

    HX = mod.(hcat(A,  B),  2)
    HZ = mod.(hcat(B', A'), 2)

    # Compute the logical operators of the code.
    (logical_X_operators, logical_Z_operators) = compute_logical_operators(HX, HZ)

    n_qubits = 2 * l * m
    k_logical = size(logical_X_operators, 1) # number of logical qubits is the number of logical operators in either X or Z basis, which should be the same.

    # Create a BBCode struct to hold the code information.
    code = BBCode(HX, HZ, logical_X_operators, logical_Z_operators, n_qubits, k_logical, distance, l, m, A_monomials, B_monomials)
    
    # Verify the commutation relations of the logical operators with the stabilizers and with each other.
    verify_canonical_commutations(HX, HZ)
    
    return code
end

function save_bb_code(bbc::BBCode; prefix::String="./../code")
    """
    Save a BB code's two matrices to plain-text files:
    <prefix>_HX.txt
    <prefix>_HZ.txt
    We will save these files in <prefix>/bb_code/, so the full paths will be `<prefix>/bb_code/<prefix>_HX.txt` and `<prefix>/bb_code/<prefix>_HZ.txt`.
    We will also write the
    - code parameters (n, k, d) in a file `<prefix>/bb_code/parameters.txt` for reference.
    - code hyperparameters (l, m, A_monomials, B_monomials) in a file `<prefix>/bb_code/hyperparameters.txt` for reference.
    - logical operators (LX, LZ) in files `<prefix>/bb_code/LX.txt` and `<prefix>/bb_code/LZ.txt` for reference.
    """
    if !isdir("$(prefix)")
        mkdir("$(prefix)")
    end
    writedlm("$(prefix)/HX.txt", bbc.HX, ' ')
    writedlm("$(prefix)/HZ.txt", bbc.HZ, ' ')
    writedlm("$(prefix)/LX.txt", bbc.LX, ' ')
    writedlm("$(prefix)/LZ.txt", bbc.LZ, ' ')

    open("$(prefix)/parameters.txt", "w") do io
        println(io, "n: ", bbc.n)
        println(io, "k: ", bbc.k)
        println(io, "d: ", bbc.d)
    end
    open("$(prefix)/hyperparameters.txt", "w") do io
        println(io, "l: ", bbc.l)
        println(io, "m: ", bbc.m)
        # Store the monomials in a human-readable format, e.g. "x^3 y^0" for the monomial (3, 0).
        println(io, "A(x, y): ", join(["x^$(a) y^$(b)" for (a, b) in bbc.A_monomials], ", "))
        println(io, "B(x, y): ", join(["x^$(a) y^$(b)" for (a, b) in bbc.B_monomials], ", "))
    end
    println("BB code saved to $(prefix)/")
    return nothing
end