# time testing
function plot_results(n_list, timings)
  figure(1)
  title("Speed of EECNetworkImpedance")  
  xlabel("n (nxn image)")
  ylabel("t [s]")
  plot(Int.(n_list), timings, "-x")
end


function elapsed_test(;n_max=70, mode_3D=false)
  timings = []
  
  n_list = collect(5 : 5 : n_max)
  for n in n_list
      println(" --- n = $(n) ---")
      push!(timings, 
          @elapsed image_to_EIS(Int64.(mode_3D ? zeros(n,n,n) : zeros(n,n)), f_list=[0.01, 100000], pyplot=false)
      )
  end
  plot_results(n_list, timings)
  return n_list, timings
end

function plot_fit(optim_list, model = (x,t) -> x[1]*exp(x[2] * t), unknown_count = 2; 
                  extrapolation_list = nothing,
                  plot_data=true,
                  plot_fit=true,
                  label=""
                  )
          
  model_zde = x -> [model(x,t) for t in optim_list[1]]
    
  function to_opt(x)
          sum = 0          
          for i in 1:length(optim_list[1])
              sum += (optim_list[2][i] - model_zde(x)[i])^2
          end
          return sum
  end
  xm = optimize(to_opt, zeros(unknown_count) .+ 1.0).minimizer
  
  if typeof(extrapolation_list) != Nothing
    plot_x_list = extrapolation_list    
  else
    plot_x_list = optim_list[1]    
  end
  plot_model = x -> [model(x,t) for t in plot_x_list]
  
  if plot_data
    plot(optim_list..., "x-")
  end
  if plot_fit
    plot(plot_x_list, plot_model(xm), label=label)
  end
  
  return plot_x_list, plot_model(xm)
end


function get_plot_for_times_fillin()
  Ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150];
  times2 = [8.618276877, 2.4021253, 1.747620653, 4.23329095, 8.971466441, 14.394831117, 
            21.057839342, 30.124806847, 45.4032383, 60.849876155, 88.453158614, 621.639396868, 
            169.945182755, 185.422826541, 241.316241434];
  times3 = [0.147667471, 2.54396585, 1.50403213, 3.478046798, 5.93846077, 10.319718565, 
            16.756460149, 25.711286162, 36.09036492, 48.51659034, 395.137235912, 99.331202957, 
            285.229349163, 780.570341478, 945.472877815];
  times4 = [3.883598193, 0.765538725, 2.426612835, 6.06551546, 13.813451025, 31.793265064, 
            44.698515326, 58.326820459, 86.772006839, 125.879504294, 167.589399414, 233.297587663, 
            293.876049016, 390.684255463, 446.523239517];

EECNetworkImpedance.plot_fit((Ns[3:15], times4[3:15]),
  #(x,t) -> x[1]*exp(x[2]*t),  2,
  #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
  #(x,t) -> x[1]*t^4, 1,
  (x,t) -> x[1]*t^3, 1,
  plot_data=true,
  label = "fill_in = 15"
  )
  # time: 2023-03-17 10:43:07 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:11], times2[3:11]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=false,
    plot_fit=true,
    extrapolation_list=Ns[3:15],
    label = "fill_in = 10"
  )
  # time: 2023-03-17 10:43:28 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:10], times3[3:10]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=false,
    plot_fit=true,
    extrapolation_list=Ns[3:15],
    label = "fill_in = 7"
  )
  # time: 2023-03-17 10:43:36 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:15], times3[3:15]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=true,
    plot_fit=false,
    #extrapolation_list=Ns[3:15],
    #label = "fill_in = 7"
  )
  # time: 2023-03-17 10:43:48 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:15], times2[3:15]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=true,
    plot_fit=false,
    #extrapolation_list=Ns[3:15],
    #label = "fill_in = 7"
  )
  # time: 2023-03-17 10:43:56 CET
  # mode: julia
  legend()
  # time: 2023-03-17 10:44:13 CET
  # mode: julia
  xlabel("N")
  # time: 2023-03-17 10:44:25 CET
  # mode: julia
  ylabel("t [s]")
  # time: 2023-03-17 10:45:05 CET
  # mode: julia
  title("study w.r.t. fill_in_ratio -> stability vs speed")
  # time: 2023-03-17 10:45:38 CET
  # mode: julia
  title("study w.r.t. fill_in_ratio with fits t=a*N^3 -> stability vs speed")

end


function MODVAL_get_plot_for_times_fillin()
  Ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150];
  times2 = [8.618276877, 2.4021253, 1.747620653, 4.23329095, 8.971466441, 14.394831117, 
            21.057839342, 30.124806847, 45.4032383, 60.849876155, 88.453158614, 621.639396868, 
            169.945182755, 185.422826541, 241.316241434];
  times3 = [0.147667471, 2.54396585, 1.50403213, 3.478046798, 5.93846077, 10.319718565, 
            16.756460149, 25.711286162, 36.09036492, 48.51659034, 395.137235912, 99.331202957, 
            285.229349163, 780.570341478, 945.472877815];
  times4 = [3.883598193, 0.765538725, 2.426612835, 6.06551546, 13.813451025, 31.793265064, 
            44.698515326, 58.326820459, 86.772006839, 125.879504294, 167.589399414, 233.297587663, 
            293.876049016, 390.684255463, 446.523239517];

  EECNetworkImpedance.plot_fit((Ns[3:15], times4[3:15]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=true,
    label = "fill_in = 15"
  )
  # time: 2023-03-17 10:43:07 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:11], times2[3:11]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=false,
    plot_fit=true,
    extrapolation_list=Ns[3:15],
    label = "fill_in = 10"
  )
  # time: 2023-03-17 10:43:28 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:10], times3[3:10]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=false,
    plot_fit=true,
    extrapolation_list=Ns[3:15],
    label = "fill_in = 7"
  )
  # time: 2023-03-17 10:43:36 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:15], times3[3:15]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=true,
    plot_fit=false,
    #extrapolation_list=Ns[3:15],
    #label = "fill_in = 7"
  )
  # time: 2023-03-17 10:43:48 CET
  # mode: julia
  EECNetworkImpedance.plot_fit((Ns[3:15], times2[3:15]),
    #(x,t) -> x[1]*exp(x[2]*t),  2,
    #(x,t) -> x[1]*t^3 + x[2]*t^2, 2,
    #(x,t) -> x[1]*t^4, 1,
    (x,t) -> x[1]*t^3, 1,
    plot_data=true,
    plot_fit=false,
    #extrapolation_list=Ns[3:15],
    #label = "fill_in = 7"
  )
  # time: 2023-03-17 10:43:56 CET
  # mode: julia
  legend()
  # time: 2023-03-17 10:44:13 CET
  # mode: julia
  xlabel("N")
  # time: 2023-03-17 10:44:25 CET
  # mode: julia
  ylabel("t [s]")
  # time: 2023-03-17 10:45:05 CET
  # mode: julia
  title("study w.r.t. fill_in_ratio -> stability vs speed")
  # time: 2023-03-17 10:45:38 CET
  # mode: julia
  title("study w.r.t. fill_in_ratio with fits t=a*N^3 -> stability vs speed")

end













function homogenous_matrix(par_study_prms::Dict)
  return generate_matrix(
            par_study_prms["dimensions"],
            par_study_prms["porosity"],
            par_study_prms["LSM_ratio"]
  )
end

function pore_shooting_matrix(par_study_prms::Dict)
  return generate_matrix(
            par_study_prms["dimensions"],
            par_study_prms["porosity"],
            par_study_prms["LSM_ratio"],
            pore_cavitance=par_study_prms["pore_cavitance"],
            LSM_cavitance=par_study_prms["LSM_cavitance"]
  )
end

function three_column_domain_matrix(p::Dict)
  return generate_matrix(
    three_column_domain_template(p["LSM_ratio1"], p["LSM_ratio2"], p["LSM_ratio3"], 
                              #
                              porosity1=p["porosity1"], porosity2=p["porosity2"], porosity3=p["porosity3"],
                              #                              
                              positions_of_contacts=p["positions_of_contacts"], height_of_contacts=p["height_of_contacts"], 
                              #
                              column_width = p["column_width"],
                              #
                              height = p["height"]
    )
  )
end

function three_column_domain_LSM_ratios(p::Dict)
  return generate_matrix(
    three_column_domain_template(p["LSM_ratios"]..., 
                              #
                              porosity1=p["porosity1"], porosity2=p["porosity2"], porosity3=p["porosity3"],
                              #                              
                              positions_of_contacts=p["positions_of_contacts"], height_of_contacts=p["height_of_contacts"], 
                              #
                              column_width = p["column_width"],
                              #
                              height = p["height"]
    )
  )
end

# function three_column_domain_LSM_ratios(p::Dict)  
#   return generate_matrix(
#     three_column_domain_template(p["LSM_ratios"]..., 
#                               #
#                               porosity1=p["porosity1"], porosity2=p["porosity2"], porosity3=p["porosity3"],
#                               #                              
#                               positions_of_contacts=p["positions_of_contacts"], height_of_contacts=p["height_of_contacts"], 
#                               #
#                               column_width = p["column_width"],
#                               #
#                               height = p["height"]
#     )
#   )
# end





















#
#
function for_each_prms_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else      
      if typeof(prms_lists[active_idx]) <: Array
        list_for_iteration = prms_lists[active_idx]
      else  
        list_for_iteration = [prms_lists[active_idx]]
      end
      for parameter in list_for_iteration
        recursive_call(push!(deepcopy(output_set),parameter), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end




function par_study(
  # if value is array, the parametric study will iterate over it  
  input_prms_dict = Dict();
  save_to_file = ""
  )
  
  function extract_physical_parameters(par_study_prms)
    parameters = []
    for physical_prms_key in filter( 
          key -> length(key) > 2 && key[1:2] == "p.",
          collect(keys(par_study_prms))
        )
      push!(parameters, physical_prms_key[3:end] => par_study_prms[physical_prms_key])
    end   
    return parameters
  end

  if typeof(input_prms_dict) <: Vector
    input_prms_dict = Dict(input_prms_dict)
  end
  
  
  par_study_prms = Dict{String, Any}(
    "repetition_idx" => 1.0,
  )
  
  if !haskey(input_prms_dict, "matrix_template")
    merge!(par_study_prms,
      Dict(
        "porosity" => 0.2,
        "LSM_ratio" => 0.2,    
        "dimensions" => (5, 5),
        "matrix_template" => homogenous_matrix
      )
    )
  end
  
  merge!(par_study_prms, input_prms_dict)

  output_data_frame = DataFrame()
  
  function execute_for_specific_values(specific_values)
    local_par_study_prms = Dict(collect(keys(par_study_prms)) .=> specific_values)     
    local_parameters = extract_physical_parameters(local_par_study_prms)
          
    if haskey(local_par_study_prms, "f_list")
      f_list = local_par_study_prms["f_list"]
    else
      f_list = "TPE"
    end

    if haskey(local_par_study_prms, "tau")
      tau = local_par_study_prms["tau"]
    else 
      tau = 1e-2
    end

    if haskey(local_par_study_prms, "verbose")
      verbose = local_par_study_prms["verbose"]
    else
      verbose = false
    end

    if haskey(local_par_study_prms, "L_el_mat")
      L_el_mat = local_par_study_prms["L_el_mat"]
    else
      L_el_mat = i_LSM
    end

    if haskey(local_par_study_prms, "R_el_mat")
      R_el_mat = local_par_study_prms["R_el_mat"]
    else
      R_el_mat = i_LSM
    end

    etime = @elapsed begin 
        mat = local_par_study_prms["matrix_template"](local_par_study_prms)
        if get(local_par_study_prms, "compute_tortuosity", false)
          @show mat
          tor = EECNetworkImpedance.matrix_tortuosity(mat)
        end
        R, R_pol, C_pol = convert.(Float64,
          image_to_EIS(
                        mat,
                        local_parameters,
                        #
                        f_list=f_list,
                        tau=tau,
                        verbose=verbose,
                        return_R_RC=true, 
                        TPE_warning=false,                      
                        pyplot=false,
                        L_el_mat=L_el_mat,
                        R_el_mat=R_el_mat,
          )
        )
      end
    
    append!(output_data_frame, 
      merge(
        local_par_study_prms,
        Dict("R" => R, "R_pol" => R_pol, "C_pol" => C_pol, "etime" => etime),
        get(local_par_study_prms, "compute_tortuosity", false) ? Dict("tor" => tor) : Dict()
      )
    )  
  end
  
  if haskey(input_prms_dict, "make_node_busy") && input_prms_dict["make_node_busy"] == 1
    dummy = 0
    while true
      #println(dummy)
      dummy += 1
    end 
  else
    for_each_prms_in_prms_lists(
      collect(values(par_study_prms)), 
      execute_for_specific_values
    )
  end
  
  #@show output_data_frame
  if save_to_file != ""
    slash_idxs = findall(x -> x in ['/'], save_to_file)
    if length(slash_idxs) > 0
      mkpath(save_to_file[1:slash_idxs[end]])
    end      
    CSV.write(save_to_file, output_data_frame)
  end
  return output_data_frame
end











































function template_par_study_three_domain()
  run_par_study(
    par_study_prms_dict = Dict(
                    "matrix_template" => EECNetworkImpedance.three_column_domain_LSM_ratios,
                #
                    "repetition_idx" => collect(1:1),
                    #
                    #
                    "LSM_ratios" => [                                     
                                    INPUT1 .* INPUT2 for INPUT1 in collect(0.0 : 0.5 : 0.0), INPUT2 in [(1, 1, 1), (0, 1, 1), (0, 0, 1)]
                                    ],
                    #
                    "porosity" => collect(0.0 : 0.05 : 0.1),
                    #
                    "positions_of_contacts" => (2, 7),
                    "height_of_contacts" => 1,
                    #
                    "column_width" => 1,
                    "height" => 10
    ),
    scripted_prms_names = [
                    "LSM_ratios", 
                    "porosity" => ["porosity1", "porosity2", "porosity3"]
    ],
    save_to_file_prefix = "3_domain_",
    direct = false,
    shell_command = "echo"
  )
end


function template_par_study_homogenous_matrix()
  EECNetworkImpedance.run_par_study(
    par_study_prms_dict = Dict(
                            "matrix_template" => EECNetworkImpedance.homogenous_matrix,
                            #
                            "repetition_idx" => collect(1:1),
                            #
                            "LSM_ratio" => collect(0.0 : 0.5 : 1.0),                            
                            "porosity" => collect(0.0 : 0.5 : 0.5),
                            #
                            "dimensions" => (5,5),                                                        
                        ), 
    scripted_prms_names = ["LSM_ratio", "porosity"],
    save_to_file_prefix = "homog_",
    direct = false,
    shell_command = "echo"
  )
end


function template_par_study_temperatue()
  for T in [600, 650, 700, 750, 800]
    for add_rep in collect(1:14)
      EECNetworkImpedance.run_par_study(
        par_study_prms_dict = Dict(
                                "matrix_template" => EECNetworkImpedance.homogenous_matrix,
                                #
                                "repetition_idx" => collect(1:3),
                                #
                                "LSM_ratio" => collect(1.0 : 0.5 : 1.0),                            
                                "porosity" => collect(0.12 : 0.001 : 0.12),
                                #
                                "dimensions" => (4,4),
                                #
                                #
                                "T" => T,
                                "p.R_YSZ" => TI("R_YSZ", T),
                                "p.R_LSM" => TI("R_LSM", T)
                            ), 
        scripted_prms_names = [],#["repetition_idx", "LSM_ratio"],
        save_to_file_prefix = "por_study_20x20_case_LSM/T$(T)_addREP$(add_rep)_",
        direct = false,
        shell_command = "echo"
      )
    end
  end
end



function template_specific_three_domain()
  # Mix ... versus ciste LSM
  mat = Dict([:M => (0.5, 0.5), :L => (0.12, 1.0)])
  
  for configuration in [(:M, :M, :M), (:L, :M, :M), (:L, :L, :M), (:L, :L, :L)]
    EECNetworkImpedance.run_par_study(
      par_study_prms_dict = Dict(
                      "matrix_template" => EECNetworkImpedance.three_column_domain_LSM_ratios,
                      # 
                      "repetition_idx" => collect(1:1),
                      #
                      #
                      "configuration" => configuration,
                      #
                      "porosity1" => mat[configuration[1]][1],
                      "porosity2" => mat[configuration[2]][1],
                      "porosity3" => mat[configuration[3]][1],
                      #
                      "LSM_ratios" => ( 
                                        mat[configuration[1]][2],
                                        mat[configuration[2]][2],
                                        mat[configuration[3]][2]                      
                                      ),                                            
                      #
                      "positions_of_contacts" => (2, 7),
                      "height_of_contacts" => 1,
                      #
                      "column_width" => 1,
                      "height" => 10,
                      #
                      #
                      "T" => T,
                      "p.R_YSZ" => TI("R_YSZ", T),
                      "p.R_LSM" => TI("R_LSM", T)
      ),
      scripted_prms_names = [
                      "LSM_ratios",                      
      ],
      save_to_file_prefix = "3_domain_",
      direct = false,
      shell_command = "echo"
    )
  end
end

function run_par_study(;shell_command="echo",
                        script_file="run_EX3_ysz_fitting_EECNetworkImpedance.jl",
                        mode="go!", 
                        direct=false,
                        par_study_prms_dict::Dict,
                        scripted_prms_names::Array=[],
                        save_to_file_prefix = "default_",
                        num_nodes=1
                      )
  
  
  
  scripted_prms_lists = [                          
                          typeof(prm_name) <: Pair ?                           
                            par_study_prms_dict[prm_name[1]] :
                            par_study_prms_dict[prm_name]                          
                          for prm_name in scripted_prms_names
                        ]
                        
  function run_par_study_for_one_dict(prms_tuple)                
    DICT = deepcopy(par_study_prms_dict)
    update_pairs = []
    for i in 1:length(scripted_prms_names)      
      if typeof(scripted_prms_names[i]) <: Pair        
        push!(update_pairs, [identifier => prms_tuple[i] for identifier in scripted_prms_names[i][2]]...)
        delete!(DICT, scripted_prms_names[i][1])
      else
        push!(update_pairs, scripted_prms_names[i] => prms_tuple[i])
      end
    end      
    DICT = merge(DICT, Dict(update_pairs))
    DICT_str = string(DICT)
    save_to_file = save_to_file_prefix*"$(prms_tuple).csv"
    #    
    if mode == "go!" || mode == "try_one_prms"
      if direct
        par_study(DICT, save_to_file=save_to_file)
      else      
        run(`$(shell_command) $(script_file) $(DICT_str) $(save_to_file) $(num_nodes)`)
      end
    elseif mode == "only_save_image"        
        @show DICT
        matrix_to_file("save_first_image.png", DICT["matrix_template"](DICT))
    end
  end
  
  if mode == "only_save_image" || mode == "try_one_prms"
    @show scripted_prms_lists
    first_prms_array = [list[Int(floor(end/2)+1)] for list in scripted_prms_lists]
    @show first_prms_array
    run_par_study_for_one_dict(first_prms_array)
  else
    for_each_prms_in_prms_lists(scripted_prms_lists, run_par_study_for_one_dict)
  end
  println("Done :) ")
  return
end


















































#### Evaluation

function test_get_processed_df()
  x_axis = "T"
  other_parameters = [
    "pore_ratio" => [0.1, 0.2, 0.3],
    "LSM_ratio"
  ]
  return show_plots(x_axis, other_parameters, "snehurka/par_studies/por_study_50x50/")  
end

function collect_df_files_in_folder(dir)
  collected_df = DataFrame()
  for file_name in readdir(dir)
    actual_DF = DataFrame(CSV.File(dir*"/$(file_name)"))
    #@show actual_DF
    #for key in ["etime", "verbose", "tau"]
    #  if key in names(actual_DF)
    #    select!(actual_DF, Not([Symbol(key)]))
    #  end
    #end
    collected_df = vcat(collected_df, actual_DF)
  end
  return collected_df
end

function test_if_all_changing_prms_were_assigned(df, keys)
  for group in groupby(df, [keys..., "repetition_idx"])
    if size(group)[1] != 1
      return false
    end
  end
  return true
end

function get_grouped_processed_df(collected_df, x_axis, other_parameters; specific_symbol="", throw_exception=true, show_var)
  processed_df = deepcopy(collected_df)    
  
  keys_to_group = []
  for other_parameter in other_parameters  
    if typeof(other_parameter) <: Pair
      filter!(row -> row[other_parameter[1]] in other_parameter[2], processed_df)
      push!(keys_to_group, other_parameter[1])
    else
      push!(keys_to_group, other_parameter)
    end
  end
  
  gdf = groupby(processed_df, [keys_to_group..., x_axis]) 
  #
  if !test_if_all_changing_prms_were_assigned(processed_df, [keys_to_group..., x_axis])
    println("\n <<< ! ! ! ERROR! not all changing parameters were specified! >>> \n")
    if throw_exception
      return throw(Exception)
    end
  end
  #
#   if show_var
    df_to_plot = combine(gdf, 
      ["R", "R_pol", "C_pol"] .=> mean, 
      ["R", "R_pol", "C_pol"] .=> var,
      ["R"] .=> maximum,
      ["R"] .=> minimum
    )
#   else
#     df_to_plot = combine(gdf, 
#       ["R", "R_pol", "C_pol"] .=> mean     
#     )
#   end
  return groupby(df_to_plot, keys_to_group)
end


function show_plots(x_axis, other_parameters, dir="snehurka/par_study/"; 
        specific_symbol="", 
        apply_func= x -> x, 
        func_name="",
        throw_exception=true,
        plot_bool=true,
        show_var=false,
        preprocessing_callback! =identity       
        )
  collected_df = collect_df_files_in_folder(dir)
  preprocessing_callback!(collected_df)
  #@show collected_df
  grouped_df = get_grouped_processed_df(
                    collected_df, 
                    x_axis, 
                    other_parameters, 
                    specific_symbol=specific_symbol, 
                    throw_exception=throw_exception, 
                    show_var=show_var
                )
  
  for (key, sub_df) in pairs(grouped_df)
    legend_entry = "$(key)"[12:end-1]
    sort!(sub_df, x_axis)
    if plot_bool
      @show sub_df
      plot_par_study_results(
        sub_df[!, x_axis], 
        sub_df[!, :R_mean], 
        sub_df[!, :R_pol_mean],
        sub_df[!, :C_pol_mean],
        label=legend_entry,
        x_axis_label=x_axis,
        title="Mean",
        apply_func=apply_func,
        func_name=func_name
      )
      if show_var
        plot_par_study_results(
          sub_df[!, x_axis], 
          sub_df[!, :R_var], 
          sub_df[!, :R_pol_var],
          sub_df[!, :C_pol_var],
          label=legend_entry,
          x_axis_label=x_axis,          
          fig_num=6,
          title="Var",
          apply_func=apply_func,
          func_name=func_name
        )
      end
    end
  end
  return grouped_df
end



function plot_par_study_results(x, R, Rp, Cp; label="", x_axis_label, fig_num=5, title="", apply_func= x -> x, func_name="")
  figure(fig_num)
  if length(title) > 0
    suptitle(title)
  end
  if func_name!=""
    ylabel_prefix = func_name*" "
  else
    ylabel_prefix = ""
  end
  
  subplot(221)
  @show R
  plot(x, apply_func.([sum(R[n])/length(R[n]) for n in 1:length(R)]), label=label, "-x")
  xlabel(x_axis_label)
  ylabel(ylabel_prefix*"R_ohm")
  legend()
  
  subplot(223)    
  plot(x, apply_func.([sum(Rp[n]/length(Rp[n])) for n in 1:length(R)]), label=label, "-x")
  xlabel(x_axis_label)
  ylabel(ylabel_prefix*"R_pol")
  legend()
  
  subplot(222)  
  plot(x, apply_func.([sum(Cp[n]/length(Cp[n])) for n in 1:length(R)]), label=label, "-x")
  xlabel(x_axis_label)
  ylabel(ylabel_prefix*"C_pol")
  legend()
    
  subplot(224)  
  plot(x, apply_func.([
      sum(R[n])/length(R[n]) + sum(Rp[n]/length(Rp[n]))
    for n in 1:length(R)]) , label=label, "-x")
  xlabel(x_axis_label)
  ylabel(ylabel_prefix*"R_tot")
  legend()
end















######  time_test_evaluation
function evaluate_time_test(dir)
  collected_df = DataFrame(N=[], times=[])
  number_of_times = 0
  for file_name in readdir(dir, join=true)
    open(file_name) do file
      N=-1
      times = []
      for ln in eachline(file)
        split_arr = split(ln, "= ")
        if occursin("eval(Meta.parse(ARGS[1]))", split_arr[1])
          prms_dict = eval(Meta.parse(split_arr[2]))
          N = prms_dict["N"]
        elseif occursin("times", split_arr[1])
          times = eval(Meta.parse(split_arr[2]))
          if number_of_times == 0
            number_of_times = length(times)-1
          elseif number_of_times > 0
            if number_of_times != length(times)-1
              return println("Error: Length of times does not match: $(number_of_times) != $(length(times))")
            end
          end
        end
      end
      if N > 0 && length(times) > 0
        #@show N, times, file_name
        push!(collected_df, (N, mean(times[2:end])))
      end
    end
  end
  sort!(collected_df)
  gdf = groupby(collected_df, "N")
  out_df = combine(gdf, :times => mean)
  return out_df
end

function plot_time_test(dir)
  df = evaluate_time_test(dir)
  plot(df[:, :N], df[:, :times_mean])
  title("Time test for NxNxN ... DIRECT LU solver")
  xlabel("N")
  ylabel("t [s]")
  return df
end

