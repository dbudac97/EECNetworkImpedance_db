module EECNetworkImpedance
using ImageIO
using FileIO
using PyPlot
using Dates
using Printf
using Colors
using DataFrames
using CSV
using Interpolations
using Optim
using SparseArrays
using Statistics
using LinearAlgebra

using IterativeSolvers
using IncompleteLU
#using ILUZero
#using Krylov

using Base

# specific material dependences

include("temperature_dependences.jl")
export TI
export TI_ID2
export TI_2D
export TI_YSZ_hi

include("physical_properties.jl")
export matrix_to_file, i_LSM, i_YSZ
export subimages_composition

# general part
include("tau_estimation.jl")
include("material_matrix_to_lin_sys.jl")
include("material_matrix_to_impedance.jl")
export material_matrix_to_impedance

# user requested part

include("equivalent_circuit_support.jl")
include("Z_view_export.jl")
include("image_to_EIS_interface.jl")
export image_to_EIS

# maybe independent package for proper geometry generation

include("generate_matrix.jl")
export generate_matrix
export three_column_domain_template
export three_column_domain_matrix
export chess_matrix
export izolated_material_ratio, characterize_material_izolation
export make_shaded_view, sym_log
export enlarge_image
export rotate_matrix, rotate_matrix!

# macro stuff consisting of par study, evaluate results, ploting, cluster operations

include("macro_stuff.jl")
export par_study


# using SymPy
# include("analytical_solution_testing.jl")
# export compare_mat
# export sym_get_A_b_from_matrix

end # module EECNetworkImpedance

