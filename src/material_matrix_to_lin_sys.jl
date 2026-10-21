# idx in A_aux matrix
# assuming that x, y can be in [1, 2, ... m*n+3]
function aux_matrix_connectivity_entries_to_lin_idx(x,y, m,n)
    if (x, y) == (1, 2)
      return 1
    elseif x == 2
      # y - 2 + 1
      return y - 1
    elseif y == m*n + 3
      #          y                - 2 + 1 ... for special case of right edge 
      # (m*n + 3) + mod(x - 3, m) - 2 + 1      
      return m*n + 2 + mod(x - 3, m)
    # general case within the matrix
    # ... horizontal ... indexed first
    elseif mod(x - y, m) == 0
      return y - 1
    # ... vertical 
    else
      # max index for horiz. case  +  idx of x  -  col_number (because there are 1 less vertical connections than #of rows)
      # (m*n + 3) + m - 2          +  x - 2     -  div(x-3, n)
      return m*(n + 1) + x - 1 - div(x-3, m)
    end
end

# 
function aux_matrix_connectivity_entries_to_lin_idx(x,y, m,n,s)
    if (x, y) == (1, 2)
      return 1
    elseif x == 2
      # left end of electrode -> horizontal
      # idx of y          + transpose to layer
      # y - 2 + 1         + m*(n+1)*div(y-3, m*n)
      return y - 1 + m*div(y-3, m*n)
    elseif y == m*n*s + 3
      # right end of electrode -> horizontal
      # 1 + #of standard horizontal of the first layer    + #row                + transpose to layer
      # 1 + m*n                                           + mod(x-3, m) + 1     + m*(n + 1)*div(x-3, m*n)
      return m*n + 2 + mod(x-3,m) + m*(n+1)*div(x-3, m*n)
    elseif div(x-3, m*n) == div(y-3, m*n) 
      if mod(x - y, m) == 0
        # horizontal
        # idx of y          + transpose to layer
        # y - 2 + 1         + m*(n+1)*div(x-3, m*n)
        return y - 1 + m*div(y-3, m*n)
      else
        # vertical
        # 1 + #horizontal    + idx of x    - #col   (because there is one connection less than row entries in each column)
        # 1 + m*(n+1)*s     + x - 2       - div(x-3, m)
        return m*(n+1)*s - 1 + x - div(x-3, m)
      end          
    else
      # interlayer
      # 1 + # horizontal  + #vertical     + idx of x
      # 1 + m*(n+1)*s     + (m-1)*n*s     + x - 2
      return (m*(2*n+1) - n)*s + x - 1     
    end
end



function ans_322(in_f, x, y)
  f(x, y) = in_f(x, y, 3,2,2)
  
  if (x,y) == (1,2); return f(x,y) == 1
  elseif (x,y) == (2,3); return f(x,y) == 2
  elseif (x,y) == (2,4); return f(x,y) == 3
  elseif (x,y) == (2,5); return f(x,y) == 4
  #
  elseif (x,y) == (3,6); return f(x,y) == 5
  elseif (x,y) == (4,7); return f(x,y) == 6
  elseif (x,y) == (5,8); return f(x,y) == 7
  #
  elseif (x,y) == (6,15); return f(x,y) == 8
  elseif (x,y) == (7,15); return f(x,y) == 9
  elseif (x,y) == (8,15); return f(x,y) == 10
  #
  #
  elseif (x,y) == (2,9); return f(x,y) == 11
  elseif (x,y) == (2,10); return f(x,y) == 12
  elseif (x,y) == (2,11); return f(x,y) == 13
  #
  elseif (x,y) == (9,12); return f(x,y) == 14
  elseif (x,y) == (10,13); return f(x,y) == 15
  elseif (x,y) == (11,14); return f(x,y) == 16
  #
  elseif (x,y) == (12,15); return f(x,y) == 17
  elseif (x,y) == (13,15); return f(x,y) == 18
  elseif (x,y) == (14,15); return f(x,y) == 19
  else
    return "ach jo"
  end
end

max_lin_idx(m,n) =  if m==1
                      aux_matrix_connectivity_entries_to_lin_idx(m*n + 2, m*n + 3, m, n)
                    else
                      aux_matrix_connectivity_entries_to_lin_idx(m*n + 1, m*n + 2, m, n)                      
                    end
max_lin_idx(m,n,s) =  if s==1
                        max_lin_idx(m,n)
                      else
                        aux_matrix_connectivity_entries_to_lin_idx(m*n*(s-1) + 2, m*n*s + 2, m, n, s)
                      end



















get_previous_list(m, n) = [(0, -1), (-1, 0)]
get_next_list(m, n) = [(1, 0), (0, 1)]

get_previous_list(m, n, s) = [(0, -1, 0), (-1, 0, 0), (0, 0, -1)]
get_next_list(m, n, s) = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]



function previous(i, j, k, auxilary_A)
  if j == 0
    return [1]
  else
    res = Array{Int64}(undef, length(get_previous_list(0,0,0)))
    for (l, d) in enumerate(get_previous_list(0,0,0))        
      res[l] = auxilary_A[(i+1, j+1, k+1) .+ d...]                      
    end      
    deleteat!(res, findall(x->x==-1,res))      
    return res
  end
end

function previous(i, j, auxilary_A)    
  if j == 0
    return [1]
  else
    res = Array{Int64}(undef, length(get_previous_list(0,0)))
    for (k, d) in enumerate(get_previous_list(0,0))        
      res[k] = auxilary_A[(i+1, j+1) .+ d...]                      
    end      
    deleteat!(res, findall(x->x==-1,res))      
    return res
  end
end

function next(i,j,k, auxilary_A)
  if j == 0          
    return auxilary_A[2:end-1, 2, 2:end-1][:]
  else
    res = Array{Int64}(undef, length(get_next_list(0,0,0)))        
    for (l, d) in enumerate(get_next_list(0,0,0))      
      res[l] = auxilary_A[(i+1, j+1, k+1) .+ d...]        
    end          
    deleteat!(res, findall(x->x==-1,res))
    return res
  end
end

function next(i,j, auxilary_A)
  if j == 0      
    return auxilary_A[2:end-1, 2]
  else
    res = Array{Int64}(undef, length(get_next_list(0,0)))
    for (k, d) in enumerate(get_next_list(0,0))
      res[k] = auxilary_A[(i+1, j+1) .+ d...]        
    end      
    deleteat!(res, findall(x->x==-1,res))
    return res
  end
end

function add_row_to_sparse_input!(sparse_input,
                                  new_row, 
                                  row_idx)        
  for item in new_row
    push!(sparse_input[1], row_idx)
    push!(sparse_input[2], item[1])
    push!(sparse_input[3], typeof(item[2]) <: Number ? float(item[2]) : item[2])
  end
end


const shift_new_row_by = 2

function add_value_to_matrix_row(matrix_header, equation_row, unknown_id, value)
    # " ..  " because U1  (in auxilary matrix) is not unknown 
    push!(matrix_header, (unknown_id-shift_new_row_by, (unknown_id)))
    push!(equation_row, (unknown_id-shift_new_row_by, value))
end

function add_equation_I_row(pos, Y_vector, auxilary_A, dims, matrix_header, sys_row_idx, sparse_input, RHS)
  new_row = []
  act_id = auxilary_A[pos .+ 1 ...]               
  #
  Y_sum = 0 + 0*im
  for node_id in previous(pos..., auxilary_A)
    Y_ji = Y_vector[aux_matrix_connectivity_entries_to_lin_idx(node_id, act_id, dims...)]
    Y_sum += Y_ji
    add_value_to_matrix_row(matrix_header, new_row, 
                            node_id,
                            Y_ji)
  end
  for node_id in next(pos..., auxilary_A)    
    Y_ij = Y_vector[aux_matrix_connectivity_entries_to_lin_idx(act_id, node_id, dims...)]
    Y_sum += Y_ij
    add_value_to_matrix_row(matrix_header, new_row, 
                            node_id, 
                            Y_ij)
  end
  add_value_to_matrix_row(matrix_header, new_row, 
    act_id,
    -Y_sum
  )

  sort!(new_row, by = first)

  # if U1 is present in the current row, it should be replaced with boundary condition U1 = 1.0
  # and removed from unknowns and pushed to RHS with opposite sign
  # ... and indexes are shifted "-1" while storing in new_row
  if new_row[1][1] == 0
    push!(RHS, -new_row[1][2])
    new_row_start_index = 2
  else
    push!(RHS, 0.0)
    new_row_start_index = 1
  end

  # if the last entry to new row is U(m*n*s)+3... the last "electrode" node... it has value 0.0 
  # ... and indexes are shifted "-1" while storing in new_row
  if new_row[end][1] == prod(dims) + 3 - shift_new_row_by
    new_row_ending_index = length(new_row) - 1
  else
    new_row_ending_index = length(new_row)
  end

  sys_row_idx[1] += 1
  add_row_to_sparse_input!(sparse_input, new_row[new_row_start_index : new_row_ending_index], sys_row_idx[1]) 
  return
end







function get_readable_matrix_header(matrix_header)
  sort!(matrix_header, by = first)
  res = [matrix_header[1]]  
  for item in matrix_header
    if item[1] != res[end][1]
      push!(res, item)
    end
  end
  return [item[2] for item in res]
end


# auxilary_A has 2 on left, m*n+3 on the right, -1 on top and bottom
function vector_to_lin_sys(Y_vector, auxilary_A)  
  # sparse_input =  (     sparse_input_row_idxs, 
  #                       sparse_input_col_idxs,
  #                       sparse_input_vals )
  matrix_header = []
  
  sparse_input = (Int64[], Int64[], ComplexF64[])
  sys_row_idx = [0]
  
  RHS = ComplexF64[]

  dims = size(auxilary_A) .- 2
  
  mode_3D = length(dims) == 3


  # other vertices
  if mode_3D
    for k in 1:dims[3], j in 1:dims[2], i in 1:dims[1]
      add_equation_I_row((i,j, k), Y_vector, auxilary_A, dims, matrix_header, sys_row_idx, sparse_input, RHS)
    end
  else
    for j in 1:dims[2], i in 1:dims[1]
        add_equation_I_row((i,j), Y_vector, auxilary_A, dims, matrix_header, sys_row_idx, sparse_input, RHS)
    end
  end
  
  matrix_header_output = get_readable_matrix_header(matrix_header)
  return matrix_header_output, sparse_input, RHS
end



































function get_auxilary_graph_matrix(m, n)
    A_aux = Array{Int64}(undef, m+2, n+2)
    #
    aux_vec = Vector{Int64}(undef, m*n)
    for i in 1:length(aux_vec)
        aux_vec[i] = i+2
    end
    pre_A = reshape(aux_vec, m, n)
    #
    A_aux[2:end-1, 2:end-1] = pre_A
    A_aux[:, 1] .= 2
    A_aux[:, end] .= m*n+3
    A_aux[1, :] .= -1
    A_aux[end, :] .= -1    
    return A_aux
end


function get_auxilary_graph_matrix(m, n, s)
    A_aux = Array{Int64}(undef, m+2, n+2, s+2)
    #
    aux_vec = Vector{Int64}(undef, m*n*s)
    for i in 1:length(aux_vec)
        aux_vec[i] = i+2
    end
    pre_A = reshape(aux_vec, m, n, s)
    #
    A_aux[2:end-1, 2:end-1, 2:end-1] = pre_A
    A_aux[:, 1, :] .= 2
    A_aux[:, end, :] .= m*n*s+3
    A_aux[1, :,  :] .= -1
    A_aux[end, :, :] .= -1
    A_aux[:, :, 1] .= -1
    A_aux[:, :, end] .= -1
    return A_aux
end

function get_large_material_matrix(pre_A::Array{<:Integer, 2}, L_el_mat, R_el_mat)
  (m,n) = size(pre_A)
  A = Array{Int64}(undef, m+2, n+2)
  A[2:end-1, 2:end-1] = pre_A
  A[:, 1] .= L_el_mat
  A[:, end] .= R_el_mat
  A[1, :] .= -1
  A[end, :] .= -1              
  return A
end

function get_large_material_matrix(pre_A::Array{<:Integer, 3}, L_el_mat, R_el_mat)            
    (m,n,s) = size(pre_A)
    A = Array{Int64}(undef, m+2, n+2, s+2)
    A[2:end-1, 2:end-1, 2:end-1] = pre_A
    A[:, 1, :] .= L_el_mat
    A[:, end, :] .= R_el_mat
    A[1, :, :] .= -1
    A[end, :, :] .= -1        
    A[:, :, 1] .= -1
    A[:, :, end] .= -1
    return A
end

function add_Y_to_vector!(Y_vector, pos, dir, aux_A, large_material_A, p, dims, w)
    n1 = large_material_A[pos...]
    n2 = large_material_A[pos .+ dir...]
    if n2 == -1
        return
    end
    
    idx1 = aux_A[pos...]
    idx2 = aux_A[pos .+ dir...]
    if idx2 < idx1
        aux = idx1
        idx1 = idx2
        idx2 = aux
    end
    
    current_key = aux_matrix_connectivity_entries_to_lin_idx(idx1, idx2, dims...)
    if Y_vector[current_key] == -Inf
        Y_vector[current_key] = 
            get_Y_entry_from_material_matrix_codes(n1, n2, p)(w)           
    else      
        Y_vector[current_key] = 1/(
          1/Y_vector[current_key] 
          +
          1/(get_Y_entry_from_material_matrix_codes(n1, n2, p)(w))
        )
    end
end










# function current_measurement(aux_A::Array{<:Integer, 2}, Y_vector, dims)
#   Y_and_act_id_sum_list = []
#   for i in 2:dims[1]+1
#     act_id = aux_A[i,2]
#     push!(Y_and_act_id_sum_list, (
#       Y_vector[aux_matrix_connectivity_entries_to_lin_idx(2, act_id, dims...)],
#       act_id  
#       )
#     )
#   end
#   return (Us, w) -> -sum([
#     (Us[act_id - shift_new_row_by] - 1.0)*Y_act for (Y_act, act_id) in Y_and_act_id_sum_list
#   ])
# end



function integrate_between_columns(col, aux_A, Y_vector, dims)
  function get_right_Us(Us, index)
    if index == 0
      return 1.0
    elseif index == length(Us)+1
      return 0.0
    else
      return Us[index]
    end
  end

  Y_and_act_id_sum_list = []
  if length(size(aux_A)) == 2
    for i in 2:dims[1]+1
      act_id = aux_A[i,col+1]
      previous_id = aux_A[i, col]
      push!(Y_and_act_id_sum_list, (
        Y_vector[aux_matrix_connectivity_entries_to_lin_idx(previous_id, act_id, dims...)],
        previous_id,
        act_id
        )
      )
    end
  else
    for i in 2:dims[1]+1, j in 2:dims[3]+1
      act_id = aux_A[i,col+1,j]
      previous_id = aux_A[i, col, j]
      push!(Y_and_act_id_sum_list, (
        Y_vector[aux_matrix_connectivity_entries_to_lin_idx(previous_id, act_id, dims...)],
        previous_id,
        act_id
        )
      )
    end
  end
  return (Us, w) -> -sum([
    (get_right_Us(Us, act_id - shift_new_row_by) - get_right_Us(Us, previous_id - shift_new_row_by))*Y_act 
    
    for (Y_act, previous_id, act_id) in Y_and_act_id_sum_list
  ])
end

function current_measurement(aux_A, Y_vector, dims)
  left = integrate_between_columns(1, aux_A, Y_vector, dims)
  right = integrate_between_columns(dims[2] + 1, aux_A, Y_vector, dims)
  return function (Us, w; side=:left)
    side == :left && return left(Us, w)
    side == :right && return right(Us, w)
    throw(ArgumentError("side must be :left or :right"))
  end
end























function material_matrix_to_lin_sys(
            material_A::Array{<:Integer}=[0 0 0; 1 1 1],
            p::parameters=parameters(),
            w::Float64=1.0;
            L_el_mat::Int64,
            R_el_mat::Int64
    )
    
    dims = size(material_A)
    mode_3D = length(dims) == 3
    aux_A = get_auxilary_graph_matrix(dims...)
    large_material_A = get_large_material_matrix(material_A, L_el_mat, R_el_mat)

    Y_vector = Vector{ComplexF64}(undef, max_lin_idx(dims...))
    Y_vector .= -Inf
    for d in vcat(get_previous_list(dims...), get_next_list(dims...))
      for i in 2:dims[1]+1
        for j in 2:dims[2]+1                
          if mode_3D
            for k in 2:dims[3]+1
              add_Y_to_vector!(Y_vector, (i, j, k), d, aux_A, large_material_A, p, dims, w)
            end
          else
            add_Y_to_vector!(Y_vector, (i, j), d, aux_A, large_material_A, p, dims, w)
          end
        end
      end
    end
    for i in 2:dims[1]+1
      if mode_3D
        for k in 2:dims[3]+1
          add_Y_to_vector!(Y_vector, (i,1, k), (0, +1, 0), aux_A, large_material_A, p, dims, w)
          add_Y_to_vector!(Y_vector, (i,dims[2]+2, k), (0, -1, 0), aux_A, large_material_A, p, dims, w)
        end
      else
        add_Y_to_vector!(Y_vector, (i,1), (0, +1), aux_A, large_material_A, p, dims, w)
        add_Y_to_vector!(Y_vector, (i,dims[2]+2), (0, -1), aux_A, large_material_A, p, dims, w)
      end
    
    end
    #return Y_vector, aux_A, large_material_A
    return vector_to_lin_sys(Y_vector, aux_A), current_measurement(aux_A, Y_vector, dims)
end
