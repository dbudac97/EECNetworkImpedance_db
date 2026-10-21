function get_prms_from_pairs(pairs)
  p = parameters()
  for pair in pairs
    setfield!(p, Symbol(pair[1]), float(pair[2]))
  end
  return p
end

function fill_static_entries!(list_eval, list)
  for i in eachindex(list)
    if !(typeof(list[i]) <: Function)
      list_eval[i] = list[i]      
    end
  end
end

function evaluate_nonstatic_entries_for_w!(eval, template, w)
  for i in eachindex(template)   
    if typeof(template[i]) <: Function
      eval[i] = template[i](w)
    end
  end
end

function set_initial_values(dims)
  x = Array{Float64}(undef, prod(dims))
  
  for k in 1:dims[3]
      for j in 1:dims[2]
          start = (k-1)*dims[1]*dims[2] + (j-1)*dims[1] + 1
          ending = start + dims[1] - 1
          x[start : ending] .= 1.0 - j/(dims[2]+1)
          #@show start, ending
      end
  end
  return x
end

function try_convert_system_to_real(sp_input, RHS)
  sp_input_values = []
  b_eval = []
  try
    sp_input_values = Float64.(sp_input[3])
    b_eval = Float64.(RHS)
  catch
    sp_input_values = sp_input[3]
    b_eval = RHS
  end
  return sp_input_values, b_eval
end

function validate_linear_solution(A, b, x; rtol=1e-8, atol=0.0)
  all(isfinite, x) || error("linear solver returned nonfinite potentials")
  residual = norm(A*x - b)
  residual <= atol + rtol*norm(b) ||
    error("linear solver residual $residual exceeds tolerance $(atol + rtol*norm(b))")
  return residual
end

function validate_electrode_currents(left, right; rtol=1e-6, atol=1e-12)
  isfinite(left) && isfinite(right) || error("nonfinite electrode current")
  abs(left - right) <= atol + rtol*max(abs(left), abs(right)) ||
    error("electrode currents do not balance: left=$left, right=$right")
  iszero(left) && error("zero total current; impedance is undefined")
  return nothing
end

function iterative_network_solution(A, b, LU; rtol, atol, max_mv_products)
  try
    x = bicgstabl(A, b, 2; Pl=LU, reltol=rtol/10,
                  max_mv_products=max_mv_products)
    validate_linear_solution(A, b, x; rtol=rtol, atol=atol)
    return x, :bicgstabl
  catch err
    # BiCGStab(l)'s small least-squares system can break down, including when
    # an accurate ILU makes the original system especially easy to solve.
    # Retry with bounded restarted GMRES, never a memory-heavy direct fallback.
    (err isa SingularException || err isa ErrorException) || rethrow()
    x = gmres(A, b; Pl=LU, reltol=rtol/10,
              restart=min(30, size(A, 1)), maxiter=max_mv_products)
    validate_linear_solution(A, b, x; rtol=rtol, atol=atol)
    return x, :gmres
  end
end

function material_matrix_to_impedance(
            material_matrix::Array = [1 1 0 1; 0 1 0 1; 0 1 0 1; 0 1 0 1],
            prms_pairs = [];
            f_list,
            iterative_solver = "auto",
            verbose = false,
            return_only_linsys = false,
            tau = "auto",
            fill_in_ratio = 10,
            solver_rtol = 1e-8,
            solver_atol = 0.0,
            current_rtol = 1e-6,
            current_atol = 1e-12,
            max_mv_products::Integer = 2000,
            solver_diagnostics = nothing,
            L_el_mat::Int64=i_LSM,
            R_el_mat::Int64=i_LSM)
  ndims(material_matrix) in (2, 3) && all(>(0), size(material_matrix)) ||
    throw(ArgumentError("expected a nonempty 2D or 3D material array"))
  all(x -> x in (i_YSZ, i_LSM, i_pore), material_matrix) ||
    throw(ArgumentError("unknown material code"))
  for value in (solver_rtol, solver_atol, current_rtol, current_atol)
    isfinite(value) && value >= 0 || throw(ArgumentError("tolerances must be finite and nonnegative"))
  end
  solver_rtol > 0 || throw(ArgumentError("solver_rtol must be positive"))
  max_mv_products > 0 || throw(ArgumentError("max_mv_products must be positive"))
  isempty(f_list) && throw(ArgumentError("frequency list must not be empty"))
  all(f -> f isa Real && !isnan(f) && f >= 0, f_list) ||
    throw(ArgumentError("frequencies must be nonnegative (0 and Inf are supported)"))
  if iterative_solver == "auto"
    iterative_solver = length(material_matrix) >= 15^3
  end
  iterative_solver isa Bool || throw(ArgumentError("iterative_solver must be true, false, or auto"))
  params = get_prms_from_pairs(prms_pairs)
  Z_list = ComplexF64[]
  for f in f_list
    verbose && @show f
    w = Float64(2*pi*f)
    (header, sp_input, RHS), measurement = material_matrix_to_lin_sys(
      material_matrix, params, w; L_el_mat=L_el_mat, R_el_mat=R_el_mat)
    values, b = try_convert_system_to_real(sp_input, RHS)
    A = sparse(sp_input[1], sp_input[2], values)
    return_only_linsys && return (A, b, measurement)
    if iterative_solver
      tau == "auto" && (tau = get_estimation_of_tau(A; target_ratio=fill_in_ratio))
      LU = ilu(A, τ=tau)
      x, method = iterative_network_solution(A, b, LU; rtol=solver_rtol,
        atol=solver_atol, max_mv_products=max_mv_products)
    else
      x = A \ b
      method = :direct
    end
    residual = validate_linear_solution(A, b, x; rtol=solver_rtol, atol=solver_atol)
    left = measurement(x, w)
    right = measurement(x, w; side=:right)
    validate_electrode_currents(left, right; rtol=current_rtol, atol=current_atol)
    impedance = 1.0/left
    isfinite(impedance) || error("nonfinite impedance at frequency $f")
    if solver_diagnostics !== nothing
      push!(solver_diagnostics, (frequency=f, method=method, residual=residual,
        left_current=left, right_current=right, current_error=abs(left-right)))
    end
    push!(Z_list, impedance)
  end
  return f_list, Z_list
end
