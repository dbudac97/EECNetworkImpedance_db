using EECNetworkImpedance
using Test

function test_matrix_3x3x3_0_0()
    PRMS = ["R_YSZ" => 10.0, "R_LSM" => 1.0, "R_pol_LSM" => 66.0, "C_pol_LSM" => 0.005]
    res3D = image_to_EIS(generate_matrix((3,3,3), 0.0, 0.0), PRMS, return_R_RC=true, TPE_warning=false, pyplot=false)
    res2D = image_to_EIS(generate_matrix((3,3), 0.0, 0.0), PRMS, return_R_RC=true, TPE_warning=false, pyplot=false)
    # the coorect answer (R_ohm, R, C) 
    the_answer = (10.33333333333333, 44.00000000000011, 0.007499999999999982)
    
    pass_3D = prod(
        isapprox.(res3D, the_answer, rtol=1e-4)
    )

    pass_2D = prod(
        isapprox.(res2D, the_answer, rtol=1e-4)
    )
    return prod([pass_2D, pass_3D])
end

function test_matrix_20x20x20_0_0()
    PRMS = ["R_YSZ" => 10.0, "R_LSM" => 1.0, "R_pol_LSM" => 66.0, "C_pol_LSM" => 0.005]
    # the coorect answer (R_ohm, R, C) 
    the_answer = (10.05, 6.6, 0.05)
    
    res2D = image_to_EIS(generate_matrix((20,20), 0.0, 0.0), PRMS, return_R_RC=true, TPE_warning=false, pyplot=false)
    res3D = image_to_EIS(generate_matrix((20,20,20), 0.0, 0.0), PRMS, return_R_RC=true, TPE_warning=false, pyplot=false)
    
    pass_3D = prod(
        isapprox.(res3D, the_answer, rtol=1e-4)
    )

    pass_2D = prod(
        isapprox.(res2D, the_answer, rtol=1e-4)
    )
    return prod([pass_2D, pass_3D])
end	

function test_chess_matrices()
    PRMS = ["R_YSZ" => 10.0, "R_LSM" => 1.0, "R_pol_LSM" => 66.0, "C_pol_LSM" => 0.005]
    # the coorect answer (R_ohm, R, C) 
    the_answer_3D = (5.446625500787345, 64.2290059945072, 0.005136950008394067)
    res3D = image_to_EIS(chess_matrix(30, 30, 30), PRMS, return_R_RC=true, TPE_warning=false, pyplot=false)
    
    the_answer_2D = (5.487020166949348, 65.5380949869996, 0.00503502392014762)
    res2D = image_to_EIS(chess_matrix(100, 100), PRMS, return_R_RC=true, TPE_warning=false, pyplot=false)

    pass_3D = prod(
        isapprox.(res3D, the_answer_3D, rtol=1e-4)
    )

    pass_2D = prod(
        isapprox.(res2D, the_answer_2D, rtol=1e-4)
    )
    return prod([pass_2D, pass_3D])
end

@testset "testing R_RC case" begin
    @test test_matrix_3x3x3_0_0()
    @test test_matrix_20x20x20_0_0()
    @test test_chess_matrices()
end


include("regressions.jl")
