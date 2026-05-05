using BicycleBivariateCodes

#============================================================================
Summary of code parameters and constructions from Table 3 of https://arxiv.org/abs/2308.07915.
Columns:
- [n, k, d] : Code parameters
- r         : Net encoding rate
- (l, m)    : Construction parameters
- A, B      : Polynomial definitions

| [n, k, d]        | r     | (l, m) | A                  | B                      |
|------------------|-------|--------|--------------------|------------------------|
| [72, 12, 6]      | 1/12  | (6, 6) | x^3 + y + y^2      | y^3 + x + x^2          |
| [90, 8, 10]      | 1/23  | (15, 3)| x^9 + y + y^2      | 1 + x^2 + x^7          |
| [108, 8, 10]     | 1/27  | (9, 6) | x^3 + y + y^2      | y^3 + x + x^2          |
| [144, 12, 12]    | 1/24  | (12, 6)| x^3 + y + y^2      | y^3 + x + x^2          |
| [288, 12, 18]    | 1/48  | (12,12)| x^3 + y^2 + y^7    | y^3 + x + x^2          |
| [360, 12, ≤24]   | 1/60  | (30, 6)| x^9 + y + y^2      | y^3 + x^25 + x^26      |
| [756, 16, ≤34]   | 1/95  | (21,18)| x^3 + y^10 + y^17  | y^5 + x^3 + x^19       |

Notes:
- Distances with ≤ indicate lower bounds.
- Polynomials A and B define the structure used in the construction.
=# # ============================================================================

function bb_code_72(;prefix::String="./../data/72q_BB_code")
    """
    Construct the [[72, 12, 6]] BB code from Table 3 of https://arxiv.org/abs/2308.07915.
    """
    bbc = gen_bb_code(6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]; distance=6)
    save_bb_code(bbc; prefix=prefix)
    return nothing
end

function bb_code_90(;prefix::String="./../data/90q_BB_code")
    """
    Construct the [[90, 8, 10]] BB code from Table 3 of https://arxiv.org/abs/2308.07915.
    """
    bbc = gen_bb_code(15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)]; distance=10)
    save_bb_code(bbc; prefix=prefix)
    return nothing
end

function bb_code_108(;prefix::String="./../data/108q_BB_code")
    """
    Construct the [[108, 8, 10]] BB code from Table 3 of https://arxiv.org/abs/2308.07915.
    """
    bbc = gen_bb_code(9, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]; distance=10)
    save_bb_code(bbc; prefix=prefix)
    return nothing
end

function bb_code_144(;prefix::String="./../data/144q_BB_code")
    """
    Construct the [[144, 12, 12]] BB code from Table 3 of https://arxiv.org/abs/2308.07915.
    """
    bbc = gen_bb_code(12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]; distance=12)
    save_bb_code(bbc; prefix=prefix)
    return nothing
end

function bb_code_288(;prefix::String="./../data/288q_BB_code")
    """
    Construct the [[288, 12, 18]] BB code from Table 3 of https://arxiv.org/abs/2308.07915.
    """
    bbc = gen_bb_code(12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)]; distance=18)
    save_bb_code(bbc; prefix=prefix)
    return nothing
end

function bb_code_360(;prefix::String="./../data/360q_BB_code")
    """
    Construct the [[360, 12, ≤24]] BB code from Table 3 of https://arxiv.org/abs/2308.07915.
    """
    bbc = gen_bb_code(30, 6, [(9, 0), (0, 1), (0, 2)], [(0, 3), (25, 0), (26, 0)]; distance=24)
    save_bb_code(bbc; prefix=prefix)
    return nothing
end

function bb_code_756(;prefix::String="./../data/756q_BB_code")
    """
    Construct the [[756, 16, ≤34]] BB code from Table 3 of https://arxiv.org/abs/2308.07915.
    """
    bbc = gen_bb_code(21, 18, [(3, 0), (0, 10), (0, 17)], [(0, 5), (3, 0), (19, 0)]; distance=34)
    save_bb_code(bbc; prefix=prefix)
    return nothing
end