using LinearAlgebra
using SparseArrays

const EEC = EECNetworkImpedance

@testset "Bounded generation and whole-packing connectivity" begin
    for cavitance in (nothing, 0.0, 0.8)
        stats = Dict{String, Any}()
        @test_throws ErrorException generate_matrix((2, 3), 1.0, 0.5;
            pore_cavitance=cavitance, max_attempts=2, generation_stats=stats)
        @test stats["attempts"] == 2
        @test stats["rejected"] == 2
        packing = generate_matrix((2, 3), 1.0, 0.5;
            pore_cavitance=cavitance, check_connectivity=false, generation_stats=stats)
        @test all(==(2), packing)
        @test stats["attempts"] == 1
        @test stats["rejected"] == 0
        @test stats["connected"] == false
    end
    for shape in ((4, 4), (4, 4, 4))
        clustered = generate_matrix(shape, .25, .5; pore_cavitance=.8,
            LSM_cavitance=.8, check_connectivity=false)
        @test count(==(2), clustered) == length(clustered) ÷ 4
        @test count(==(1), clustered) == 3length(clustered) ÷ 8
    end
    # An intentionally pore-only block is subsequently overwritten by contacts.
    template = [[(1, 1), (3, 1), 1.0, 0.5],
                [(1, 2), (3, 4), 0.0, 0.5],
                [(2, 1), (2, 1), 0.0, 1.0]]
    packing = generate_matrix(template; max_attempts=1)
    @test EEC.check_material_connection(packing)
    @test packing[1, 1] == packing[3, 1] == 2
    @test packing[2, 1] == 1
    disconnected = [[(1, 1), (3, 4), 1.0, 0.5]]
    @test_throws ErrorException generate_matrix(disconnected; max_attempts=2)
    @test all(==(2), generate_matrix(disconnected; check_connectivity=false))
    @test_throws ArgumentError generate_matrix((2, 3), -0.1, 0.5)
    @test_throws ArgumentError generate_matrix((2, 3), 0.0, 0.5; max_attempts=0)
    @test size(generate_matrix((2, 3), 0.0, 0.5; recursion_depth=1)) == (2, 3)
    @test_throws ErrorException generate_matrix((2, 3), 1.0, 0.5; recursion_depth=2)
end

@testset "Rotations preserve orientation and support rectangular domains" begin
    square = [1 2; 3 4]
    original = copy(square)
    @test rotate_matrix(square, 1, "z") == [3 1; 4 2]
    @test square == original
    @test rotate_matrix!(square, 1, "z") === square
    @test square == [3 1; 4 2]
    rectangular = [1 2 3; 4 5 6]
    @test rotate_matrix(rectangular, 1, "z") == [4 1; 5 2; 6 3]
    @test_throws DimensionMismatch rotate_matrix!(rectangular, 1, "z")
    @test rectangular == [1 2 3; 4 5 6]
    @test rotate_matrix!(copy(rectangular), 2, "z") == [6 5 4; 3 2 1]
    cube = reshape(collect(1:8), 2, 2, 2)
    for (axis, plane) in (("x", (2, 3)), ("y", (1, 3)), ("z", (1, 2)))
        expected = similar(cube)
        for idx in CartesianIndices(cube)
            src = collect(Tuple(idx))
            src[plane[1]] = 3 - idx[plane[2]]
            src[plane[2]] = idx[plane[1]]
            expected[idx] = cube[src...]
        end
        @test rotate_matrix(cube, 1, axis) == expected
        @test rotate_matrix!(copy(cube), 1, axis) == expected
    end
    box = reshape(collect(1:24), 2, 3, 4)
    for (axis, shape) in (("x", (2, 4, 3)), ("y", (4, 3, 2)), ("z", (3, 2, 4)))
        rotated = rotate_matrix(box, 1, axis)
        @test size(rotated) == shape
        @test rotate_matrix(rotated, -1, axis) == box
        @test rotate_matrix(box, 4, axis) == box
        @test_throws DimensionMismatch rotate_matrix!(copy(box), 1, axis)
        @test rotate_matrix!(copy(box), 2, axis) == rotate_matrix(box, 2, axis)
    end
end

@testset "Solid shortest paths and isolation" begin
    @test EEC.matrix_tortuosity(fill(0, 2, 5)) == 1.0
    @test EEC.matrix_tortuosity(fill(0, 2, 5, 3)) == 1.0
    # The left pore must not be treated as a valid origin of a solid path.
    @test EEC.matrix_tortuosity([2 0 0]) == Inf
    @test EEC.matrix_tortuosity(fill(2, 2, 3)) == Inf
    detour = [0 2 0 0; 0 0 0 2]
    @test EEC.matrix_tortuosity(detour) > 1
    @test characterize_material_izolation([0 0 0; 2 2 2]) == [4 4 4; 0 0 0]
    @test characterize_material_izolation([0 2 1]) == [2 0 3]
    @test characterize_material_izolation([2 0 2]) == [0 1 0]
    @test izolated_material_ratio([2 0 2]) == 1.0
end

@testset "ILU search terminates for unattainable fill" begin
    @test EEC.get_estimation_of_tau(spdiagm(0 => ones(10)); target_ratio=12) == 0.0
    @test EEC.get_estimation_of_tau(sparse(reshape([2.0], 1, 1)); target_ratio=12) == 0.0
    A = spdiagm(-1 => -ones(9), 0 => fill(3.0, 10), 1 => -ones(9))
    @test isfinite(EEC.get_estimation_of_tau(A; target_ratio=0.1, max_evaluations=4))
    @test_throws ArgumentError EEC.get_estimation_of_tau(A; target_ratio=0)
    # Larger samples must search positive tolerances rather than factor at tau=0.
    @test EEC.get_estimation_of_tau(spdiagm(0 => ones(400));
        target_ratio=12, max_evaluations=4) > 0
end

@testset "Residual and electrode-current validation rejects bad solutions" begin
    A = spdiagm(0 => ones(2))
    @test EEC.validate_linear_solution(A, [1., 2.], [1., 2.]) == 0
    @test_throws ErrorException EEC.validate_linear_solution(A, [1., 2.], [0., 0.])
    @test_throws ErrorException EEC.validate_linear_solution(A, [1., 2.], [NaN, 0.])
    @test_throws ErrorException EEC.validate_electrode_currents(1., 2.)
    @test_throws ErrorException EEC.validate_electrode_currents(NaN, 1.)
end

@testset "Analytical full spectrum, electrode materials, and exact endpoints" begin
    @test EEC.RC_el(40., .005, Inf) == 0
    @test EEC.RC_el(40., 0., Inf) == 40
    @test EEC.RC_el(40., .005, 0.) == 40
    prms = ["R_YSZ" => 10., "R_LSM" => 1., "R_pol_LSM" => 66., "C_pol_LSM" => .005]
    frequencies = [0., .01, 1., 100., Inf]
    # Independent analytical result for uniform YSZ with LSM contacts.
    for shape in ((3, 5), (3, 5, 2)), iterative in (false, true)
        area = shape[1] * (length(shape) == 3 ? shape[3] : 1)
        expected = [(10*shape[2] + 1 + (isinf(f) ? 0 : 132/(1 + 2pi*im*f*66*.005)))/area
                    for f in frequencies]
        diagnostics = []
        _, actual = material_matrix_to_impedance(fill(0, shape), prms;
            f_list=frequencies, iterative_solver=iterative, tau=0., solver_diagnostics=diagnostics)
        @test actual ≈ expected rtol=1e-7
        @test length(diagnostics) == length(frequencies)
        _, ionic = material_matrix_to_impedance(fill(0, shape), prms;
            f_list=frequencies, iterative_solver=iterative, tau=0., L_el_mat=i_YSZ, R_el_mat=i_YSZ)
        @test ionic ≈ fill(10*(shape[2]+1)/area, length(frequencies)) rtol=1e-7
    end
    # Keep endpoint extraction and the original finite-frequency TPE defaults.
    r0, rp, c = image_to_EIS(fill(0, 3, 5), prms; f_list=[0., Inf], return_R_RC=true, pyplot=false)
    @test r0 ≈ 10.2
    @test rp ≈ 26.4
    @test c == 0
    @test all(isapprox.(image_to_EIS(fill(0, 3, 5), prms; f_list=[Inf, 0.], return_R_RC=true, pyplot=false), (r0, rp, c)))
    fs, zs = image_to_EIS(fill(0, 3, 5), prms; TPE_warning=false, pyplot=false)
    @test length(fs) == length(zs) == 21
    @test all(isfinite, zs)
    @test_throws ArgumentError material_matrix_to_impedance(fill(0, 2, 2); f_list=[-1.])
end

@testset "Heterogeneous full spectra agree across solvers" begin
    prms = ["R_YSZ" => 10., "R_LSM" => 1., "R_pol_YSZ" => 8.,
            "C_pol_YSZ" => .002, "R_pol_LSM" => 66., "C_pol_LSM" => .005]
    frequencies = [0., .01, 1., 100., Inf]
    packing = chess_matrix(3, 4, 2)
    _, direct = material_matrix_to_impedance(packing, prms;
        f_list=frequencies, iterative_solver=false)
    _, iterative = material_matrix_to_impedance(packing, prms;
        f_list=frequencies, iterative_solver=true, tau=0.)
    @test direct ≈ iterative rtol=1e-7
    _, rotated = material_matrix_to_impedance(rotate_matrix(packing, 1, "y"), prms;
        f_list=frequencies, iterative_solver=false)
    @test direct ≈ rotated rtol=1e-7
end
