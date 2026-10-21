
# All connectivity checks treat LSM and YSZ as a single solid phase.
function validate_generation_inputs(dimensions, porosity, LSM_ratio,
                                    pore_cavitance, LSM_cavitance)
  length(dimensions) in (2, 3) || throw(ArgumentError("expected 2D or 3D dimensions"))
  all(d -> d isa Integer && d > 0, dimensions) ||
    throw(ArgumentError("dimensions must be positive integers"))
  for (name, value) in (("porosity", porosity), ("LSM_ratio", LSM_ratio),
                        ("pore_cavitance", pore_cavitance), ("LSM_cavitance", LSM_cavitance))
    value === nothing && continue
    isfinite(value) && 0 <= value <= 1 ||
      throw(ArgumentError("$name must be between 0 and 1"))
  end
  pore_cavitance === nothing && LSM_cavitance != 0.0 &&
    throw(ArgumentError("set pore_cavitance when using LSM_cavitance"))
end

function generation_attempt_limit(max_attempts, recursion_depth)
  limit = recursion_depth === nothing ? max_attempts : recursion_depth
  isfinite(limit) && isinteger(limit) && 1 <= limit <= typemax(Int) ||
    throw(ArgumentError("max_attempts (or legacy recursion_depth) must be a positive integer"))
  return Int(limit)
end

function generate_connected_domain(make_domain; check_connectivity, max_attempts,
                                    generation_stats=nothing)
  for attempt in 1:max_attempts
    domain = make_domain()
    connected = (check_connectivity || generation_stats !== nothing) ?
                  check_material_connection(domain) : nothing
    if generation_stats !== nothing
      generation_stats["attempts"] = attempt
      generation_stats["rejected"] = check_connectivity ? attempt - (connected ? 1 : 0) : 0
      generation_stats["connected"] = connected
    end
    (!check_connectivity || connected) && return domain
  end
  throw(ErrorException("No solid path between electrodes after $max_attempts attempts; " *
                       "use check_connectivity=false to retain unconnected packings."))
end

function generate_matrix(dimensions, porosity::Float64, LSM_ratio::Float64;
    pore_cavitance::Union{Float64, Nothing}=nothing,
    LSM_cavitance::Union{Float64, Nothing}=0.0,
    recursion_depth=nothing, max_attempts::Integer=1000,
    check_connectivity=true, generation_stats=nothing)
  validate_generation_inputs(dimensions, porosity, LSM_ratio, pore_cavitance, LSM_cavitance)
  limit = generation_attempt_limit(max_attempts, recursion_depth)
  if pore_cavitance !== nothing
    return shoot_pores(dimensions, porosity, LSM_ratio;
      pore_cavitance=pore_cavitance, LSM_cavitance=LSM_cavitance,
      max_attempts=limit, check_connectivity=check_connectivity,
      generation_stats=generation_stats)
  end
  function make_domain()
    domain = Array{Int16}(undef, dimensions)
    for i in eachindex(domain)
      domain[i] = rand() < porosity ? i_pore : (rand() < LSM_ratio ? i_LSM : i_YSZ)
    end
    return domain
  end
  return generate_connected_domain(make_domain; check_connectivity=check_connectivity,
    max_attempts=limit, generation_stats=generation_stats)
end


function generate_submatrix_to_matrix(matrix, left_upper::Union{Tuple, Array}, right_lower::Union{Tuple, Array}, porosity::Float64, LSM_ratio::Float64)
  submatrix = generate_matrix(right_lower .- left_upper .+ (1,1), porosity, LSM_ratio; check_connectivity=false)
  matrix[left_upper[1] : right_lower[1], left_upper[2] : right_lower[2]] = deepcopy(submatrix)
  return 
end

function generate_submatrix_to_matrix(matrix, right_lower::Union{Tuple, Array}, porosity::Float64, LSM_ratio::Float64)
  generate_submatrix_to_matrix(matrix, (1,1), right_lower::Union{Tuple, Array}, porosity::Float64, LSM_ratio::Float64)
end

function get_default_submatrix_list()
  # 1 item in the resulting list is
  # 
  # left upper, right lower corner, porosity, LSM_ratio

  first_block_column = 1
  second_block_column = 21
  third_block_column = 41
  width = 42
  height = 80
  return [
    [(1,1), (height, first_block_column), 1.0, 0.5],
    [(1,first_block_column + 1), (height, second_block_column), 0.2, 0.5],
    [(1,second_block_column + 1), (height, third_block_column), 0.2, 1.0],
    [(1,third_block_column + 1), (height, width), 0.2, 1.0],
    #
    [(10, 1), (15, first_block_column), 0.0, 1.0],
    [(40, 1), (45, first_block_column), 0.0, 1.0],
    [(60, 1), (65, first_block_column), 0.0, 1.0]
  ]
end

function generate_matrix(submatrix_list::Array=get_default_submatrix_list();
                         check_connectivity=true, max_attempts::Integer=1000,
                         recursion_depth=nothing, generation_stats=nothing)
  limit = generation_attempt_limit(max_attempts, recursion_depth)
  return generate_connected_domain(; check_connectivity=check_connectivity,
      max_attempts=limit, generation_stats=generation_stats) do
    assemble_submatrices(submatrix_list)
  end
end

function assemble_submatrices(submatrix_list)
  isempty(submatrix_list) && throw(ArgumentError("submatrix list must not be empty"))
  max_height, max_width = -1, -1
  for submatrix in submatrix_list
    s = submatrix[2]
    if s[1] > max_height
      max_height = s[1]
    end
    if s[2] > max_width
      max_width = s[2]
    end    
  end  
  
  output_matrix = Matrix{Int16}(undef, max_height, max_width)
  output_matrix .= -1

  for submatrix in submatrix_list
    generate_submatrix_to_matrix(output_matrix, submatrix...)    
  end 
  invalid_submatrix_list = false
  foreach(x -> x==-1 ?  invalid_submatrix_list = true : false, output_matrix)
  if invalid_submatrix_list
    throw(ArgumentError("submatrices must cover every pixel"))
  else
    return output_matrix
  end
end



function three_column_domain_template(LSM_ratio1, LSM_ratio2, LSM_ratio3; 
                              #
                              porosity1=0.5, porosity2=0.5, porosity3=0.5,
                              #                              
                              positions_of_contacts=[15, 50], height_of_contacts=5, 
                              #
                              column_width = 5,
                              #
                              height = 70
                              )
  # 1 item in the resulting list is
  # 
  # left upper, right lower corner, porosity, LSM_ratio
  first_block_column = 1
  second_block_column = column_width + 1
  third_block_column = 2*column_width + 1
  fourth_block_column = 3*column_width + 1
  width = fourth_block_column + 1
  
  output = [
    [(1,1), (height, first_block_column), 1.0, 0.5],
    [(1,first_block_column + 1), (height, second_block_column), porosity1, LSM_ratio1],
    [(1,second_block_column + 1), (height, third_block_column), porosity2, LSM_ratio2],
    [(1,third_block_column + 1), (height, fourth_block_column), porosity3, LSM_ratio3],    
    [(1, fourth_block_column + 1), (height, width), 0.0, 1.0]
  ]
  
  # contacts
  for p_of_contact in positions_of_contacts
    push!(output, 
        [(p_of_contact, 1), (p_of_contact + height_of_contacts - 1, first_block_column), 0.0, 1.0]
    )
  end
  
  return output
end
function chess_matrix(m,n,s)
  A = Array{Int64}(undef, (m,n,s))
  A .= 0
  for i in 1:m, j in 1:n, k in 1:s
    if mod(i+j+k, 2) == 0
      A[i,j,k] = 1
    end
  end
  return A
end

function chess_matrix(m, n)
  A = Matrix{Int64}(undef, m,n)
  A .= 0
  for i in 1:m, j in 1:n
    if mod(i+j, 2) == 0
      A[i,j] = 1
    end
  end
  return A
end

function aux_domain(domain::Array{T} where T<: Integer; inner_number=nothing, boundary_number=0)
  if length(size(domain)) == 2
    aux_matrix = Matrix{Integer}(undef, size(domain) .+ (2,2))
    if typeof(inner_number) == Nothing
      aux_matrix[2:end-1, 2:end-1] = domain
    else
      aux_matrix[:, :] .= inner_number
    end
    aux_matrix[:, 1] .= boundary_number
    aux_matrix[:, end] .= boundary_number
    #
    aux_matrix[1, :] .= boundary_number
    aux_matrix[end, :] .= boundary_number
    #
    return aux_matrix
  elseif length(size(domain)) == 3
    aux_matrix = Array{Integer}(undef, size(domain) .+ (2,2,2))
    if typeof(inner_number) == Nothing
      aux_matrix[2:end-1, 2:end-1, 2:end-1] = domain
    else
      aux_matrix[:, :, :] .= inner_number
    end

    aux_matrix[:, end, :] .= boundary_number
    aux_matrix[:, 1, :]   .= boundary_number
    #
    aux_matrix[1, :, :]   .= boundary_number
    aux_matrix[end, :, :] .= boundary_number
    #
    aux_matrix[:, :, 1]   .= boundary_number
    aux_matrix[:, :, end] .= boundary_number
    #
    return aux_matrix
  else
    prinltn("ERROR: length(size(domain)) $(length(size(domain))) != 2 or 3")
    return throw(Exception)
  end  
end

function search_dirs(domain::Array{T} where T<: Integer)
  if length(size(domain)) == 2
    return search_dirs = [(1,0), (-1,0), (0, 1), (0, -1)]
  elseif length(size(domain)) == 3
    return search_dirs = [(1,0,0), (-1,0,0), (0, 1,0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
  else
    prinltn("ERROR: length(size(domain)) $(length(size(domain))) != 2 or 3")
    return throw(Exception)
  end
end

function prepare_aux_matrix(domain)
  if length(size(domain)) == 2
    dims = 2
  else
    dims = 3
  end  
  
  if dims == 2
    aux_matrix = Matrix{Integer}(undef, size(domain) .+ (2,2))
    aux_matrix .= 1 
    aux_matrix[:, end] .= 0
    aux_matrix[:, 1] .= 0
    #
    aux_matrix[1, :] .= 0
    aux_matrix[end, :] .= 0
    
    search_dirs = [(1,0), (-1,0), (0, 1), (0, -1)]
  end
  
  if dims == 3
    aux_matrix = Array{Integer}(undef, size(domain) .+ (2,2,2))
  
    aux_matrix .= 1 
    aux_matrix[:, end, :] .= 0
    
    aux_matrix[:, 1, :] .= 0
    #
    aux_matrix[1, :, :] .= 0
    aux_matrix[end, :, :] .= 0
    #
    aux_matrix[:, :, 1] .= 0
    aux_matrix[:, :, end] .= 0
    
    search_dirs = [(1,0,0), (-1,0,0), (0, 1,0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
  end
  return aux_matrix
end

function start_point(col, dims, gen_row)
  res = []
  push!(res, gen_row[1] + 1)
  push!(res, col)
  if dims == 3
      push!(res, gen_row[2] + 1)
  end
    
  return tuple(res...)
end


function spread_number!(aux_matrix, domain, num; from_side = "left")
  if length(size(domain)) == 2
    dims = 2
    generalized_rows = [(row, layer) for row in 1:size(domain)[1], layer in 1:1]
  else
    dims = 3
    generalized_rows = [(row, layer) for row in 1:size(domain)[1], layer in 1:size(domain)[3]]
  end
  start_col = (from_side == "left" ? 2 : size(domain)[2]+1)
  
  for generalized_row in generalized_rows
    list_to_process = [start_point(start_col, dims, generalized_row)]
    
    while length(list_to_process) > 0
      aux_coors = list_to_process[end]
      deleteat!(list_to_process, length(list_to_process))
      
      if (aux_matrix[aux_coors...] == 1)
        if (domain[(aux_coors .- 1)...] in i_material_list)          
          aux_matrix[aux_coors...] = num         
          for dir in search_dirs(domain)
            push!(list_to_process, aux_coors .+ dir)
          end
        else
          aux_matrix[aux_coors...] = 0
        end
      end
    end
  end
end

function characterize_material_izolation(domain::Array{T} where T <: Integer)
  # 0 = pore
  # 1 = izolated material
  # 2 = connected to left only
  # 3 = connected to right only
  # 4 = connected to both
  
  L_connected_sign = 4 
  L_matrix = prepare_aux_matrix(domain)
  spread_number!(L_matrix, domain, L_connected_sign, from_side="left")

  R_connected_sign = 7
  R_matrix = prepare_aux_matrix(domain)
  spread_number!(R_matrix, domain, R_connected_sign, from_side="right")

  sum_matrix = R_matrix .+ L_matrix

  res_matrix = prepare_aux_matrix(domain)
  dims = length(size(domain))
  if dims == 2
    res_matrix[2:end-1, 2:end-1] .= domain
  else
    res_matrix[2:end-1, 2:end-1, 2:end-1] .= domain
  end

  for i in 1:length(res_matrix[:])
    pix = 0
    if res_matrix[i] in i_material_list
      if sum_matrix[i] <= 2
        pix = 1
      elseif sum_matrix[i] in [4,5]
        pix = 2
      elseif sum_matrix[i] in [7,8]
        pix = 3
      elseif sum_matrix[i] in [11]
        pix = 4
      end
    end
    res_matrix[i] = pix
    
  end
 
  if dims == 2
    return res_matrix[2:end-1, 2:end-1]
  else
    return res_matrix[2:end-1, 2:end-1, 2:end-1]
  end
end

function test_izolated_material_ratio()
  mmm = [0 0 0 2; 2 2 0 2; 2 0 0 2; 2 2 2 2;;; 2 2 0 2; 2 2 0 0; 2 2 0 2; 2 2 2 2;;; 2 2 0 2; 0 0 2 0; 2 0 2 2; 2 2 2 2;;; 2 2 2 2; 2 2 0 2; 2 2 0 2; 2 2 0 2]
  return (
    izolated_material_ratio(mmm, wanted_ids = [1]) == 0.16666666666666666 &&
    izolated_material_ratio(mmm, wanted_ids = [4]) == 0.6666666666666666 &&
    izolated_material_ratio(mmm, wanted_ids = [1,2,3]) == 0.3333333333333333
  )
end

function izolated_material_ratio(domain::Array{T} where T<: Integer; wanted_ids=[1])
  char_matrix = characterize_material_izolation(domain)
  all_mat = 0
  wanted_mat = 0
  for i in 1:length(char_matrix[:])
    if char_matrix[i] > 0
      all_mat += 1
      if char_matrix[i] in wanted_ids
        wanted_mat += 1
      end
    end
  end
  return wanted_mat/all_mat
end

function check_material_connection(domain::Array{T} where T <: Integer)
  material_sign = 8
  end_sign = 6
  
  if length(size(domain)) == 2
    dims = 2
  else
    dims = 3
  end  
  
  if dims == 2
    aux_matrix = Matrix{Integer}(undef, size(domain) .+ (2,2))
    aux_matrix .= 1 
    aux_matrix[:, end] .= end_sign
    aux_matrix[:, 1] .= 0
    #
    aux_matrix[1, :] .= 0
    aux_matrix[end, :] .= 0
    
    search_dirs = [(1,0), (-1,0), (0, 1), (0, -1)]
  end
  
  if dims == 3
    aux_matrix = Array{Integer}(undef, size(domain) .+ (2,2,2))
  
    aux_matrix .= 1 
    aux_matrix[:, end, :] .= end_sign
    
    aux_matrix[:, 1, :] .= 0
    #
    aux_matrix[1, :, :] .= 0
    aux_matrix[end, :, :] .= 0
    #
    aux_matrix[:, :, 1] .= 0
    aux_matrix[:, :, end] .= 0
    
    search_dirs = [(1,0,0), (-1,0,0), (0, 1,0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
  end
  
  is_connected=false  
  
  
  if dims==2    
    for row in 1:size(domain)[1]
      list_to_process = [(row +1, 2)]
      
      while length(list_to_process) > 0
        x, y = list_to_process[end]
        deleteat!(list_to_process, length(list_to_process))
        
        if aux_matrix[x,y] == end_sign        
          is_connected=true
          return true
        else
          if (aux_matrix[x,y] == 1)
            #@show x,y
            if (domain[x-1, y-1] in i_material_list)          
              aux_matrix[x,y] = material_sign         
              for dir in search_dirs
                push!(list_to_process, (x,y) .+ dir)
                #search_around_this(((x,y) .+ dir)...)
              end
            else
              aux_matrix[x,y] = 0
            end
          end
        end
      end            
    end
  
  
  else

    
    for row in 1:size(domain)[1], layer in 1:size(domain)[3]
      list_to_process = [(row +1, 2, layer+1)]
      
      while length(list_to_process) > 0
        x, y, z = list_to_process[end]
        deleteat!(list_to_process, length(list_to_process))
        
        if aux_matrix[x,y,z] == end_sign        
          is_connected=true
          return true
        else
          if (aux_matrix[x,y,z] == 1)           
            if (domain[x-1, y-1, z-1] in i_material_list)                      
              aux_matrix[x,y,z] = material_sign         
              for dir in search_dirs        
                push!(list_to_process, (x,y,z) .+ dir)                
              end
            else
              aux_matrix[x,y, z] = 0
            end        
          end
        end
      end
      
    end    
  end
  return is_connected
end



































function typical_use_shoot_pores()
  mm = EECNetworkImpedance.shoot_pores((20, 20, 20), 0.7, 0.4, pore_cavitance=0.0, LSM_cavitance = 0.1)
  EECNetworkImpedance.matrix_to_file("jojojjoooojoooo.png", mm)
  println(EECNetworkImpedance.check_material_connection( 
      mm
      )
  )

end




function test_shoot_pores()
  return shoot_pores((10, 10), 0.5, 0.4, pore_cavitance=0.1, LSM_cavitance = 0.2)
end


function ext_correction(domain)
  if length(size(domain)) == 2
    return (1, 1)
  else
    return (1,1,1)
  end  
end

function check_item_is_boundary(ext_domain, pos, target_i)
  is_boundary = false
  ext_corr = ext_correction(ext_domain)
  if ext_domain[pos .+ ext_corr...] == -1
    return false
  end
  for dir in search_dirs(ext_domain)
    if (ext_domain[pos .+ dir .+ ext_corr...] != target_i) && 
        (ext_domain[pos .+ dir .+ ext_corr...] != i_pore) && 
        (ext_domain[pos .+ dir .+ ext_corr...] != -1)
      is_boundary = true
    end
  end
  return is_boundary
end

function get_indicies_of_all_neighbours(ext_domain, pos, target_i)
  neigh_list = []
  for dir in search_dirs(ext_domain)
    test_pos = pos .+ dir .+ ext_correction(ext_domain)
    if (ext_domain[test_pos...] != target_i) && (ext_domain[test_pos...] != -1) && (ext_domain[test_pos...] != i_pore)
      push!(neigh_list, test_pos)
    end
  end
  return neigh_list
end

function get_indicies_of_random_neighbour(ext_domain, pos, target_i)
  neigh_list = get_indicies_of_all_neighbours(ext_domain, pos, target_i)
  return neigh_list[Int32(rand(1:length(neigh_list)))] .- ext_correction(ext_domain)  
end

function get_standard_domain(ext_domain)
  if length(size(ext_domain)) == 2
    return ext_domain[2:end-1, 2:end-1]
  else
    return ext_domain[2:end-1, 2:end-1, 2:end-1]
  end
end

function get_body_list(domain)
  dims = size(domain)
  if length(dims) == 2
    return [(x,y) for x in 1:dims[1], y in 1:dims[2]][:]
  elseif length(size(domain)) == 3
    return [(x,y,z) for x in 1:dims[1], y in 1:dims[2], z in 1:dims[3]][:]
  else
    prinltn("ERROR: length(size(domain)) $(length(size(domain))) != 2 or 3")
    return throw(Exception)
  end  
end

function shoot_pores(dims, porosity, LSM_ratio; pore_cavitance, LSM_cavitance,
                     recursion_depth=nothing, max_attempts::Integer=1000,
                     check_connectivity=true, generation_stats=nothing)
  validate_generation_inputs(dims, porosity, LSM_ratio, pore_cavitance, LSM_cavitance)
  limit = generation_attempt_limit(max_attempts, recursion_depth)
  return generate_connected_domain(; check_connectivity=check_connectivity,
      max_attempts=limit, generation_stats=generation_stats) do
    shoot_pores_once(dims, porosity, LSM_ratio; pore_cavitance=pore_cavitance,
                     LSM_cavitance=LSM_cavitance)
  end
end

function shoot_pores_once(dims, porosity, LSM_ratio; pore_cavitance, LSM_cavitance)
  
  domain = Array{Integer}(undef, dims)
  body_list = get_body_list(domain)
  #
  pix_tot = prod(dims)
  extended_domain = aux_domain(domain, inner_number=i_YSZ, boundary_number=-1)
  
  # pore shooting ###########################
  pix_por = Int32(round(porosity*pix_tot))
  boundary_pore_list = []
  for i in 1:pix_por        
    mother_item_idx = nothing
    if (rand() <= pore_cavitance) && (length(boundary_pore_list) > 0)        
        mother_item_idx = rand(1:length(boundary_pore_list))        
        swapping_item_indices = get_indicies_of_random_neighbour(extended_domain, boundary_pore_list[mother_item_idx], i_pore)        
        swapping_item_idx = findall(x -> x == swapping_item_indices, body_list)[1]
    else
        swapping_item_idx = rand(1:length(body_list))
        swapping_item_indices = body_list[swapping_item_idx]
    end
    
    
    extended_domain[swapping_item_indices .+ ext_correction(domain)...] = i_pore
    
    if (typeof(mother_item_idx) != Nothing) && !check_item_is_boundary(extended_domain, boundary_pore_list[mother_item_idx], i_pore)
      deleteat!(boundary_pore_list, mother_item_idx)      
    end
    
    # delete boundary_pore_list items which happens to be interior all of the sudden[[
    for dir in search_dirs(domain)
        
        if !check_item_is_boundary(extended_domain, swapping_item_indices .+ dir, i_pore)          
          search_result = findall(x -> x == swapping_item_indices .+ dir, boundary_pore_list)
          
          if length(search_result) > 0
            deleteat!(boundary_pore_list, search_result[1])
          end
        end
    end
    
    if check_item_is_boundary(extended_domain, swapping_item_indices, i_pore)
      if !(swapping_item_indices in boundary_pore_list)
        push!(boundary_pore_list, swapping_item_indices)
      end
    end    
    deleteat!(body_list, swapping_item_idx)
  end

  #@show body_list
  # LSM shooting ###########################
  if LSM_cavitance > 0.0
    pix_LSM = Int32(round((1 - porosity)*pix_tot*LSM_ratio))
    boundary_LSM_list = []
    for i in 1:pix_LSM     
      
      mother_item_idx = nothing
      if (rand() <= LSM_cavitance) && (length(boundary_LSM_list) > 0)        
          
          mother_item_idx = rand(1:length(boundary_LSM_list))        
          swapping_item_indices = get_indicies_of_random_neighbour(extended_domain, boundary_LSM_list[mother_item_idx], i_LSM)        
          swapping_item_idx = findall(x -> x == swapping_item_indices, body_list)[1]
      else
          swapping_item_idx = rand(1:length(body_list))
          swapping_item_indices = body_list[swapping_item_idx]
      end
      
      
      extended_domain[swapping_item_indices .+ ext_correction(domain)...] = i_LSM
      
      if (typeof(mother_item_idx) != Nothing) && !check_item_is_boundary(extended_domain, boundary_LSM_list[mother_item_idx], i_LSM)
        deleteat!(boundary_LSM_list, mother_item_idx)      
      end
      
      # delete boundary_LSM_list items which happens to be interior all of the sudden[[
      for dir in search_dirs(domain)
          
          if !check_item_is_boundary(extended_domain, swapping_item_indices .+ dir, i_LSM)          
            search_result = findall(x -> x == swapping_item_indices .+ dir, boundary_LSM_list)
            
            if length(search_result) > 0
              deleteat!(boundary_LSM_list, search_result[1])
            end
          end
      end
      
      if check_item_is_boundary(extended_domain, swapping_item_indices, i_LSM)
        if !(swapping_item_indices in boundary_LSM_list)
          push!(boundary_LSM_list, swapping_item_indices)
        end
      end    
      deleteat!(body_list, swapping_item_idx)
    end

    for i in 1:length(extended_domain[:])
      if (extended_domain[i] != -1) && (extended_domain[i] != i_pore) && (extended_domain[i] != i_LSM)
        extended_domain[i] = i_YSZ
      end
    end
  else
    
    # the rest
    for i in 1:length(extended_domain[:])
      if (extended_domain[i] != -1) && (extended_domain[i] != i_pore)
        if rand() < LSM_ratio
          extended_domain[i] = i_LSM
        else
          extended_domain[i] = i_YSZ
        end
      end
    end
  end



  return get_standard_domain(extended_domain)
end








my_fu(x,a) = log(9, x*a+1);
norm_f(x, a) = my_fu(x,a)/(2*my_fu(0.5,a))
sym_log(x, a) = if x <= 0.5
                return norm_f(x,a)
            else
                return 1-norm_f(1 -x,a)
            end
#my_range = collect(0 : 0.01 : 1.0)
#plot(my_range, [sym_f(x,100) for x in my_range])







function save_palette(p1, p2, file_name)
  save(file_name, hcat(p1, p2))
end

function make_shaded_view(domain::Array{<:Integer, 3}, file_name; grad_depth=15, only_palette = false,
                            LSM_s = (0.3, 0.3, 0.3), LSM_e = (0.0, 0.0, 0.0),
                            YSZ_s = (1.0, 0.8, 0.0), YSZ_e = (0.5, 0.4, 0.0),
                            pall_func = x->x,
                            side=1
                          )

  if side in [1,2,3,4]
    axis = "x"
    step = side-1
  elseif side in [5,6]
    axis = "y"
    step = (side == 5 ? 1 : -1)
  end
  domain = rotate_matrix(domain, step, axis)

  LSM_palette = [RGB(LSM_s .+ (LSM_e .- LSM_s) .*pall_func(x) ...) for x in range(0.0, 1.0, grad_depth)]
  YSZ_palette = [RGB(YSZ_s .+ (YSZ_e .- YSZ_s) .*pall_func(x) ...) for x in range(0.0, 1.0, grad_depth)]

  if only_palette
    save_palette(LSM_palette, YSZ_palette, file_name)
    return
  end

  res_img = Array{Union{RGB, Nothing}}(undef, size(domain[:, :, 1]))
  res_img .= nothing
  for layer in 1:size(domain)[3]
    if layer > grad_depth
      cor_layer = grad_depth
    else
      cor_layer = layer
    end
    for x in 1:size(domain)[1], y in 1:size(domain)[2]
      if typeof(res_img[x,y]) == Nothing
        if domain[x,y,layer] == i_LSM
          res_img[x,y] = LSM_palette[cor_layer]
        elseif domain[x,y,layer] == i_YSZ
          res_img[x,y] = YSZ_palette[cor_layer]  
        end
      end
    end
  end

  res_img = map(x -> (typeof(x) == Nothing ? RGB(1., 1., 1.) : x), res_img)

  save(file_name, res_img)
  return 
end




function rotate_2D(M, step=1)
  return rotr90(M, mod(step, 4))
end

function rotation_axes(A, axis)
  if ndims(A) == 2
    axis in ("x", "y", "z") || throw(ArgumentError("unknown rotation axis $axis"))
    return (1, 2)
  end
  ndims(A) == 3 || throw(ArgumentError("expected a 2D or 3D array"))
  axis == "x" && return (2, 3)
  axis == "y" && return (1, 3)
  axis == "z" && return (1, 2)
  throw(ArgumentError("unknown rotation axis $axis"))
end

function rotate_matrix(A, step::Integer, axis)
  first_axis, second_axis = rotation_axes(A, axis)
  result = copy(A)
  # The same orientation as the original square-slice rotations.
  for _ in 1:mod(step, 4)
    permutation = collect(1:ndims(A))
    permutation[first_axis], permutation[second_axis] = second_axis, first_axis
    result = reverse(permutedims(result, permutation); dims=second_axis)
  end
  return result
end

function rotate_matrix!(A, step::Integer, axis)
  first_axis, second_axis = rotation_axes(A, axis)
  if isodd(step) && size(A, first_axis) != size(A, second_axis)
    throw(DimensionMismatch("a quarter-turn changes array shape; use rotate_matrix"))
  end
  copyto!(A, rotate_matrix(A, step, axis))
  return A
end


function search_wide(domain; starting_point, feasible_idxs=[], blocking_idxs=[])
  boundary_number = -1
  inner_number = typemax(Int)

  aux_domain = EECNetworkImpedance.aux_domain(domain, inner_number=inner_number, boundary_number=boundary_number)
  # podivam se na cislo ve fronte a pridam sousedni souradnice a uz jim rovnou napisu cisla
  
  starting_idx = 1
  list_to_process = [starting_point .+ 1]
  
  ending_aux_coord = nothing
  start = Tuple(starting_point)
  checkbounds(Bool, domain, start...) || throw(ArgumentError("start is outside the domain"))
  if !(domain[start...] in i_material_list)
    return aux_domain, 0
  end

  aux_domain[(starting_point .+ 1)...] = 1
  while starting_idx <= length(list_to_process)
    aux_coors = list_to_process[starting_idx]
    act_num = aux_domain[(aux_coors)...]
  
    
         
    for dir in search_dirs(domain)
      new_aux_coors = aux_coors .+ dir
      if new_aux_coors[2] == size(aux_domain)[2]
        ending_aux_coord = aux_coors
        break
      end

      if (aux_domain[new_aux_coors...] != boundary_number) &&
            (aux_domain[new_aux_coors...] > act_num +1) &&
            (domain[(new_aux_coors .- 1)...] in i_material_list)
        push!(list_to_process, new_aux_coors)
        aux_domain[new_aux_coors...] = act_num + 1
      end     
    end

    if typeof(ending_aux_coord) != Nothing
      break
    end
    starting_idx += 1    
  end
  
  
  #@show ending_aux_coord, aux_domain[ending_aux_coord...]
  if typeof(ending_aux_coord) == Nothing
    return aux_domain, 0
  else
    return aux_domain, aux_domain[ending_aux_coord...]
  end
end

function matrix_tortuosity(domain)
  paths = []
  
  for x in 1:size(domain)[1], z in 1:(  length(size(domain)) == 3 ? size(domain)[3] : 1)
      ad, p = EECNetworkImpedance.search_wide(domain; starting_point = length(size(domain)) == 3 ? [x,1,z] : [x,1])
      p > 0 ? push!(paths, p) : Nothing
  end
  # Inf explicitly denotes no solid path; do not silently average an empty set.
  return isempty(paths) ? Inf : mean(paths)/size(domain, 2)
end


