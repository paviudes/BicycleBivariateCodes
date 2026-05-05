using BicycleBivariateCodes
using DelimitedFiles

function test_commutations()
    println("Running logical operators tests...")
    
    # Load the stabilizer matrices from ./../data/tests/test_Steane.
    stabilizers_X = readdlm("./../data/tests/test_Steane/HX.txt", ' ', Int)
    stabilizers_Z = readdlm("./../data/tests/test_Steane/HZ.txt", ' ', Int)
    result_Steane = verify_canonical_commutations(stabilizers_X, stabilizers_Z)

    # Load the stabilizer matrices from ./../data/tests/test_Hamming_hgp.
    stabilizers_X = readdlm("./../data/tests/test_Hamming_hgp/HX.txt", ' ', Int)
    stabilizers_Z = readdlm("./../data/tests/test_Hamming_hgp/HZ.txt", ' ', Int)
    result_Hamming_hgp = verify_canonical_commutations(stabilizers_X, stabilizers_Z)
    
    if result_Steane && result_Hamming_hgp
        println("\n✓ All tests passed!")
        return true
    else
        println("\n✗ Some tests failed!")
        return false
    end
end