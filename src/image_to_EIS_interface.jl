
function get_material_matrix(matrix_input)
  input_path = ""
  if typeof(matrix_input) == String
    input_path = matrix_input
    material_matrix = file_to_matrix(matrix_input)
  else
    material_matrix = matrix_input
  end
  if size(material_matrix)[1] < 1 || size(material_matrix)[2] < 1
    println("ERROR: wrong input matrix")
    return
  end
  return input_path, material_matrix
end

function change_extension_to(path, new_ext)
  last_dot_position = findall(x -> x == '.', path)[end]
  if last_dot_position < 2
    println("ERROR: input_path: last_dot_position = $(last_dot_position) < 2")
  end
  return path[1:last_dot_position-1]*"."*new_ext
end

function image_to_EIS(
            matrix_input::Union{Array, String} = [1 1 0 1; 0 1 0 1; 0 1 0 1; 0 1 0 1],            
            prms_pairs = [];
            #            
            f_list="TPE",
            TPE_f_list_in = [1e-3, 1e6],
            TPE_f_list_out = [10.0^n for n in (-3 : 0.5 : 7)],
            TPE_warning = true,
            #f_list=[2.0^i for i in -3:10], 
            #
            iterative_solver = "auto",
            tau="auto",
            fill_in_ratio = 12,
            solver_rtol = 1e-8,
            solver_atol = 0.0,
            current_rtol = 1e-6,
            current_atol = 1e-12,
            max_mv_products::Integer = 2000,
            solver_diagnostics = nothing,
            verbose = false,
            pyplot=true,
            return_R_RC=false,
            export_z_file="",
            #export_z_file="!use_file_name",
            save_also_image="",
            #save_also_image="!input",
            
            
            L_el_mat=i_LSM,
            R_el_mat=i_LSM,
            
            
            return_specific_impedance = true,
            store_R_RC=""
            )
  
  input_path, material_matrix = get_material_matrix(matrix_input)
  
  two_point_extrapolation = false
  if typeof(f_list)==String && (f_list == "two_point_extrapolation" || f_list == "TPE")
    two_point_extrapolation = true
    f_list = TPE_f_list_in
  end
  
  extract_R_RC = false  
  if return_R_RC || two_point_extrapolation || store_R_RC != ""
    extract_R_RC = true
  end
                    
  f_list, Z_list =  material_matrix_to_impedance(
    material_matrix,         
    prms_pairs,
    #      
    f_list = f_list,            
    #
    iterative_solver = iterative_solver,
    verbose = verbose,
    tau = tau,
    fill_in_ratio = fill_in_ratio,
    solver_rtol=solver_rtol, solver_atol=solver_atol,
    current_rtol=current_rtol, current_atol=current_atol,
    max_mv_products=max_mv_products, solver_diagnostics=solver_diagnostics,
    #
    L_el_mat=L_el_mat,
    R_el_mat=R_el_mat
  )
  
  if return_specific_impedance
    #normalize output
    dims = size(material_matrix)
    if length(dims) == 2
      Z_list .*= dims[1] / dims[2]
    else
      Z_list .*= dims[1]*dims[3] / dims[2]
    end
  end
  
  
  if extract_R_RC 
    if length(f_list) == 1 && f_list[1] == Inf
      R_ohm, R, C = (real(Z_list[1]), 0.0, 0.0)
    elseif length(f_list) == 2 && f_list[1] == Inf && f_list[2] == 0.0
      R_ohm, R, C = (real(Z_list[1]), real(Z_list[2])-real(Z_list[1]), 0.0)
    elseif length(f_list) == 2 && f_list[1] == 0.0 && f_list[2] == Inf
      R_ohm, R, C = (real(Z_list[2]), real(Z_list[1])-real(Z_list[2]), 0.0)
    else
      R_ohm, R, C = (1.0, 1.0, 1.0)
      if length(f_list) < 2
        println("ERROR: extract R_ohm, R, C: length(f_list) < 2, $(length(f_list)) < 2")
        return
      else
        omegas = 2*pi .* [f_list[1], f_list[end]]
        Zs = [Z_list[1], Z_list[end]]
        try
          R_ohm, R, C = get_R_RC_prms_from_2Z(omegas, Zs)      
        catch        
          println("ERROR: R_RC circuit evaluation failed for (Z1, Z2) = $(Zs)  ... returning (f_list, Z_list)\n")
          extract_R_RC = false
        end      
      end
    end
  end
  
  if two_point_extrapolation && extract_R_RC    
    f_list = TPE_f_list_out
    Z_list = get_Z_from_R_RC(f_list, R_ohm, R, C)
    TPE_warning && println("WARNING: TPE (two_point_extrapolation) has been used!")
  end

  pyplot && nyquistPlot(Z_list)
  
  if store_R_RC != "" && extract_R_RC
    open(store_R_RC, "a") do f
      write(f, @sprintf("%s\t%s\t%E\t%E\t%E\n", 
        Dates.now(), 
        (input_path == "" ? "<matrix_input>" : input_path ), 
        R_ohm, R, C)
      )
    end
  end
  
  
  if export_z_file == "!use_file_name"
    if input_path != ""      
      export_EIS_to_Z_file(change_extension_to(input_path, "z"), f_list, Z_list)
    else
      println("ERROR: cannot export_z_file: input_path == \"\"!")
      return throw(Exception)
    end
    if save_also_image == "!asZfile"
      println("ERROR: cannot save image, it already exists!")
      return throw(Exception)
    end
  elseif export_z_file != ""
    
    export_EIS_to_Z_file(export_z_file, f_list, Z_list)

    if save_also_image == "!asZfile"
      if input_path !="" 
        aux_extension = input_path[end-2 : end]
      else
        aux_extension = "png"
      end
      matrix_to_file(
        change_extension_to(export_z_file, aux_extension),
        material_matrix
      )    
    end
  end
  
  if export_z_file == ""
    if save_also_image == "!asZfile" 
      println("ERROR: cannot save image sith Z file name, no Z file name defined!")
      throw(Exception)
    end
  end
  
  if (save_also_image != "")
    if (save_also_image[1] != '!')   
      matrix_to_file(save_also_image, material_matrix)
    else
      println("WARNING: be careful, you are trying to save image with a wrong name!")
    end
  end
  
  
  if extract_R_RC && return_R_RC
    return R_ohm, R, C
  else
    return f_list, Z_list  
  end
end

