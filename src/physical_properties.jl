const i_YSZ = 0
const i_LSM = 1
const i_pore = 2

const i_material_list = [0,1]


Base.@kwdef mutable struct parameters        
    R_YSZ::Float64 = 1/0.045 # S/cm
    R_pol_YSZ::Float64 = 0
    C_pol_YSZ::Float64 = 0.001
    #
    R_LSM::Float64 = 1/290 # S/cm
    R_pol_LSM::Float64 = 40
    C_pol_LSM::Float64 = 0.005
    #
    R_pore::Float64 = 1e12
end

function RC_el(R, C, w)
  # Evaluate exact endpoints explicitly, including a purely resistive C=0 branch.
  (iszero(R) || (isinf(w) && C > 0)) && return 0.0
  (iszero(C) || iszero(w)) && return R
  return R/(1 + C*R*w*im)
end

function get_Y_entry_from_material_matrix_codes(n1, n2, p::parameters)
    if      n1 == i_LSM
            if      n2 == i_LSM
              return w -> 1/(p.R_LSM/2)
            elseif  n2 == i_YSZ
              if p.R_pol_LSM == 0.0
                return w -> 1/(p.R_LSM/2)
              else
                return w -> 1/(p.R_LSM/2 + RC_el(p.R_pol_LSM, p.C_pol_LSM, w))
              end
            elseif  n2 == i_pore
              return w ->  1/(p.R_LSM/2)
            end
    elseif  n1 == i_YSZ
            if      n2 == i_LSM
              if p.R_pol_YSZ == 0.0
                return w -> 1/(p.R_YSZ/2)
              else          
                return w -> 1/(p.R_YSZ/2 + RC_el(p.R_pol_YSZ, p.C_pol_YSZ, w))
              end             
            elseif  n2 == i_YSZ
              return w -> 1/(p.R_YSZ/2)
            elseif  n2 == i_pore
              return w -> 1/(p.R_YSZ/2)
            end
    elseif  n1 == i_pore            
            return w -> 1/(p.R_pore/2)
    else
        println("ERROR: get_Y_entry...")
    end        
end



function file_to_matrix(path="src/geometry.png")
  RGB_m = load(path)
  m = Matrix{Int}(undef, size(RGB_m)...)
  for i in 1:size(m)[1]
    for j in 1:size(m)[2]
      (r,g,b) = (RGB_m[i,j].r, RGB_m[i,j].g, RGB_m[i,j].b)
      if (r > 0.5) && (g > 0.5) && (b > 0.5)
        m[i,j] = i_pore
      elseif (r > 0.5) && (g > 0.5) 
        m[i,j] = i_YSZ
      else
        m[i,j] = i_LSM
      end
    end
  end
  return m
end

function matrix_to_file(path, matrix, colors=((1.0,0.8,0), (0,0,0)))
  save(path, map(x -> if     x == i_pore
                   RGB(1, 1, 1)
                elseif x == i_YSZ
                   RGB(colors[1]...)
                elseif x == i_LSM
                   RGB(colors[2]...)
                end, matrix)
  )
  return
end

function enlarge_image(path; resolution=nothing, resize_factor=nothing, print_bool=false)
  RGB_m = load(path)
  n = size(RGB_m)[1]
  if typeof(resize_factor)==Nothing
    if typeof(resolution)==Nothing
      println("ERROR: No resize factor or resolution defined!")
      throw(Exception)
      return
    elseif typeof(resolution) <: Integer
      resize_factor = Integer(ceil(resolution/size(RGB_m)[1]))
    else
      println("ERROR: Bad type of resolution parameter!")
      throw(Exception)
    end
  end
  big_image = Array{RGB}(undef, size(RGB_m) .* resize_factor)
  for i in 1:n, j in 1:n
    big_image[
      resize_factor*(i-1) + 1 : resize_factor*(i), 
      resize_factor*(j-1) + 1 : resize_factor*(j)
    ] .= RGB_m[i,j]
  end
  save("$(path[1:end-4])_resized.$(path[end-2 : end])", big_image)
  print_bool && println("done: $(path[1:end-4])_resized.$(path[end-2 : end])")
  return
end

function name_from_template(template="first#(D)_second#(B)", dict=Dict("D"=>5, "B"=>6))
  replace_array = ["#($key)" => "$(value)" for (key, value) in dict]
  replace(template, replace_array...)
end

function subimages_composition(;folder_path, 
                                row_numbers, col_numbers, template, 
                                background_color=(1,1,1), del_width=10, out_filename="compilation.png")
  function file_xy(dict)
    specific_path = name_from_template(template, dict)
    load("$(folder_path)/$(specific_path)")
  end
  
  M = length(row_numbers)
  N = length(col_numbers)
  
  first_file = file_xy(Dict("row"=>row_numbers[1], "col"=>col_numbers[1]))
  m,n = size(first_file)

  big_image = Array{RGB}(undef, (M*m + (M-1)*del_width, N*n + (N-1)*del_width))
  big_image[1:end, 1:end] .= RGB(background_color...)

  for (x_idx, x) in enumerate(row_numbers), (y_idx, y) in enumerate(col_numbers)
    big_image[
      (m + del_width)*(x_idx - 1) + 1 : (m + del_width)*(x_idx - 1) + m,
      (n + del_width)*(y_idx - 1) + 1 : (n + del_width)*(y_idx - 1) + n,
    ] .= file_xy(Dict("row"=>x, "col"=>y))
  end
  save("$(folder_path)/$(out_filename)", big_image)
  return
end



