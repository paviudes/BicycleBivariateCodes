module BicycleBivariateCodes

using Nemo
using LinearAlgebra
using DelimitedFiles

include("commutations.jl")
export verify_canonical_commutations

include("bbcodes.jl")
export BBCode, save_bb_code, gen_bb_code

include("logicals.jl")
export nullspace_mod2, compute_logical_operators, canonicalize_logicals

end # module BicycleBivariateCodes