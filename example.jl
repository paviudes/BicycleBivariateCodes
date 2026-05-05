#!/usr/bin/env julia

"""
Example usage of the BicycleBivariateCodes.jl package.

This script demonstrates how to:
1. Load the package 
2. Generate some BB codes
3. Compute logical operators
"""

using Pkg

# Activate the local environment (this will load the Project.toml and create Manifest.toml if needed)
Pkg.activate(".")

# Load the package
using BicycleBivariateCodes

function main()
    println("=== BicycleBivariateCodes.jl Example ===")
    
    # Generate a BB code
    println("\n1. Generating BB code [72, 12, 6]...")
    code_72 = bb_code_72(prefix="./data/codes/example_72q_BB_code")
    println("Generated BB code with parameters: n=$(code_72.n), k=$(code_72.k), d=$(code_72.d)")
    
    # Compute logical operators for the generated code
    println("\n2. Computing logical operators...")
    logical_X, logical_Z = compute_logical_operators(code_72.HX, code_72.HZ)
    println("Logical X operators shape: $(size(logical_X))")
    println("Logical Z operators shape: $(size(logical_Z))")
    
    # Generate another code
    println("\n3. Generating BB code [90, 8, 10]...")
    code_90 = bb_code_90(prefix="./data/codes/example_90q_BB_code")
    println("Generated BB code with parameters: n=$(code_90.n), k=$(code_90.k), d=$(code_90.d)")
    
    println("\n=== Example completed successfully! ===")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end