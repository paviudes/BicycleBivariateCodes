# BicycleBivariateCodes

A Julia package for generating **Bivariate Bicycle (BB) codes** — quantum low-density parity-check (LDPC) codes whose distance and encoding rate both grow with the number of physical qubits.

We follow the construction in [Bravyi et al. (2024), *High-threshold and low-overhead fault-tolerant quantum memory*](https://www.nature.com/articles/s41586-024-07107-7) to build the X- and Z-type stabilizers of these CSS codes, and use Gaussian elimination over GF(2) to compute a canonical set of logical operators $\bar{X}_i$ and $\bar{Z}_i$ for $1 \leq i \leq k$.

All codes in Table 6 of the [supplementary material](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07107-7/MediaObjects/41586_2024_7107_MOESM1_ESM.pdf) and in Table 3 of [arXiv:2308.07915](https://arxiv.org/abs/2308.07915) can be generated directly with this package.

---

## Setup (first time only)

1. [Install Julia](https://julialang.org/downloads/) if you haven't already.
2. Clone the repository:
   ```
   git clone git@github.com:paviudes/BicycleBivariateCodes.git
   ```
3. Enter the repo:
   ```
   cd BicycleBivariateCodes
   ```
4. Start the Julia REPL, then activate and instantiate the package environment:
   ```julia
   julia> ]activate .
   (BicycleBivariateCodes) pkg> instantiate
   (BicycleBivariateCodes) pkg> resolve
   ```
5. Press Backspace to leave the package mode, then `exit()` to quit the REPL.

---

## Generating codes from Table 3 of arXiv:2308.07915

1. Change into the `expts/` directory and start Julia:
   ```
   cd expts/
   julia
   ```
2. Activate the package and import it:
   ```julia
   julia> ]activate ./../
   julia> using BicycleBivariateCodes
   ```
3. Call the appropriate code-generation function, supplying an output directory via `prefix`:
   ```julia
   julia> bb_code_72(; prefix="./../data/72q_BB_code")
   ```
   Replace `72` with any `n` value from Table 3 (e.g. `bb_code_90`, `bb_code_108`, `bb_code_144`, `bb_code_288`, `bb_code_360`, `bb_code_756`).

The output directory will contain:

| File | Contents |
|---|---|
| `HX.txt` | X-type stabilizers (binary matrix, one row per stabilizer) |
| `HZ.txt` | Z-type stabilizers |
| `LX.txt` | X-type logical operators |
| `LZ.txt` | Z-type logical operators |
| `parameters.txt` | Code parameters `n`, `k`, `d` |
| `hyperparameters.txt` | Generation parameters: `l`, `m`, polynomials `A(x,y)` and `B(x,y)` |

---

## API reference

### Generating a code from scratch

`gen_bb_code(l, m, A_monomials, B_monomials; distance)` constructs a BB code from the ring parameters `l`, `m` and the monomial lists of the polynomials `A(x,y)` and `B(x,y)` from Section 4 of [the paper](https://arxiv.org/abs/2308.07915). It returns a `BBCode` struct containing the stabilizers, logical operators, and code parameters.

```julia
# [[72, 12, 6]] BB code
bbc = gen_bb_code(6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]; distance=6)
```

### Saving a code to disk

`save_bb_code(bbc; prefix)` writes the stabilizers, logical operators, and parameter files to the given directory:

```julia
save_bb_code(bbc; prefix="./../data/72q_BB_code")
```

### Computing logical operators for any CSS code

`compute_logical_operators(HX, HZ)` accepts the binary parity-check matrices of any CSS code (rows = stabilizers) and returns a canonical pair `(LX, LZ)` satisfying the symplectic commutation relations — see [Wilde (2009)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.79.062322) for the algorithm:

```julia
LX, LZ = compute_logical_operators(HX, HZ)
```

### Verifying canonical commutation relations

`verify_canonical_commutations(HX, HZ, LX, LZ)` checks the canonical commutation relations (Section 3.1 of [Gottesman (2013)](https://arxiv.org/abs/1310.3235)):
- every X-stabilizer commutes with every logical Z operator,
- every Z-stabilizer commutes with every logical X operator,
- $\bar{X}_i$ anticommutes with $\bar{Z}_i$ and commutes with $\bar{Z}_j$ for $j \neq i$.

```julia
verify_canonical_commutations(HX, HZ, LX, LZ)
```
