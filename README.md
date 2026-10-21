# EECNetworkImpedance.jl

[![DOI](https://zenodo.org/badge/620190626.svg)](https://zenodo.org/badge/latestdoi/620190626)

## Introduction

The goal is to simulate an impedance measurement of an electrochemical cell of defined structure consisting of YSZ, LSM and pores using electrical elements. 

### Geometry
The structure is defined by a material matrix which contains integers

- *0* : for a YSZ particle (yellow),
- *1* : for a LSM particle (black),
- *2* : for a pore (white).

Material matrix can be also defined by a raster image using colors specified above in brackes.

!["Specified cell geometry"](images/geometry.png?raw=true )

The left and right contacts are selected independently using `L_el_mat` and `R_el_mat`. Either contact may be LSM or YSZ. YSZ contacts impose the potential through the ion-conducting material, corresponding to an idealized ionic potentiostat in this network model.

!["Cell geometry with electordes"](images/geometry_with_electrodes.png?raw=true )

### Electrochemistry
Physical simulation is done using standard electrical elements - R elements for ohmic resistance and RC element for a double-layer behavior on interfaces of LSM | YSZ. We have the following parameters:

- *R_YSZ* : ohmic resistance of YSZ
- *C_pol_YSZ* : polarization capacitance on YSZ side
- *R_pol_YSZ* : polarization resistance on YSZ side
- *R_LSM* : ohmic resistance of LSM
- *C_pol_LSM* : polarization capacitance on LSM side
- *R_pol_LSM* : polarization resistance on LSM side
- *R_pore* : resistance of pore (is set to very high number as default)


Geometry specific behavior is simulated via interactions between pixels by the manner described in the following picture.

!["Particle impedance scheme"](images/scheme.png?raw=true )


The each impedance *Z* is specified using the information about starting material *M1* and the ending material *M2* (see the picture) in the following manner:

!["Interaction scheme"](images/scheme_interaction.png?raw=true )

- *M1 = YSZ*
  - *M2 = YSZ*  => *Z = R_YSZ/2* 
  - *M2 = LSM*  => *Z = R_YSZ/2 + Z_RC(R_pol_YSZ, C_pol_YSZ)* (including the double-layer on YSZ side)
  - *M2 = pore* => *Z = R_YSZ/2*
- *M1 = LSM*
  - *M2 = YSZ*  => *Z = R_LSM/2 + Z_RC(R_pol_LSM, C_pol_LSM)* (including the double-layer on LSM side)
  - *M2 = LSM*  => *Z = R_LSM/2*
  - *M2 = pore* => *Z = R_LSM/2*
- *M1 = pore* => *Z = R_pore/2*

The factor 1/2 is there so that the total resistance through the whole (e. g.) YSZ particle is *R_YSZ*.

## Installation
The package can be then installed via 
```julialang
]
add https://github.com/dbudac97/EECNetworkImpedance_db
```


## Usage

Before using the package, you have to execute

```julialang
using EECNetworkImpedance
```

### Basics

Supposing we have either 

- `material_matrix = [1 1 1; 0 1 2]`(with values in {0,1,2}) or 
- bitmap image *my_image.png* (with colors {yellow, black, white}). 

In addition, we can specify parameters as a set of pairs. The following are the default parameters:

```julialang
physical_parameters = [ "R_YSZ" => 1/0.045,
                        "R_pol_YSZ" => 0,
                        "C_pol_YSZ" => 0.001, 
                        #
                        "R_LSM" => 1/290,
                        "R_pol_LSM" => 40, 
                        "C_pol_LSM" => 0.005, 
                        #
                        "R_pore" => 1e12]
```

TPE is a fast estimate of effective conductivity and aggregate polarization resistance for large statistical studies, not an exact reconstruction of a general network spectrum. It is particularly useful when simulating thousands of random packings. Highly organized geometries can exhibit multiple prominent arcs even with one local interfacial RC type. To study spectral shape or relaxation modes, supply an explicit frequency list and solve the full network at every requested frequency.

The default TPE inputs are finite frequencies `[1e-3, 1e6]` Hz. Exact `0.0` and `Inf` frequencies are also supported. With the pair `[0.0, Inf]` (in either order), resistance extraction returns the high-frequency resistance and the difference between the endpoint resistances; capacitance is returned as zero because these two limits do not determine a relaxation time. Use this pair for resistance extraction (`return_R_RC=true`), rather than reconstruction of a dispersive spectrum.

If less parameters are specified, the others are supposed to be default, i.e.

```julialang
physical_parameters = ["R_YSZ" => 73]
```

The core function is

```julialang
f_list, Z_list = image_to_EIS(material_matrix, physical_parameters)
```

or using path to image file

```julialang
f_list, Z_list = image_to_EIS("images/geometry.png", physical_parameters)
```

or specifying parameters in a function call 

```julialang
f_list, Z_list = image_to_EIS("images/geometry.png", ["R_YSZ" => 73])
```

or without specifying parameters


```julialang
f_list, Z_list = image_to_EIS("images/geometry.png")
```

which returns (by default) frequencies `f_list` for which impedances `Z_list` are computed.

### Additinal options

Practically useful keyword parameters are

- `f_list = [1, 10, 100]` : specification of array of frequencies for which EIS simulation will run. 
  - `= [10.0^n for n in (-3 : 0.5 : 7)]` is a good format for a geometric sequence with base 10.0 (it must be floating point number, not integer 10) and exponent spanning from -3 to 7 with a step 0.5.
  - `= "two_point_extrapolation"` : the simulation is run only for `TPE_f_list_in` yielding two impedances, 
      R-RC circuit is analytically fitted to the two computed impedances. The output Z_list is computed using this R-RC circuit for 
      frequencies in `TPE_f_list_out`.
  - default value is `= "TPE"` : which is a shortcut for "two_point_extrapolation" with the same meaning
- `TPE_f_list_in = [1e-3, 1e6]` by default
- `TPE_f_list_out = [10.0^n for n in (-3 : 0.5 : 7)]` by default
- `TPE_warning = true` : if true and `"TPE" = true`, a warning about performing TPE is printed to standard output. 
- `pyplot = true` : if *false*, no Nyquist plot is plotted
- `return_R_RC = false` :
  - if `= true` : the output of function `image_to_EIS` is a tripple (R_ohm, R_pol, C_pol) from R-RC circuit
  - if `= false` : the output is a tuple `(f_list, Z_list)`
- `export_z_file = ""` : decides whether a standard file for Z_view is exported
  - default value is `= ""`, which means *do nothing*
  - if `= "some_file.z"` : exports to this file
  - if `= "!use_file_name"` : this option is valid only when the function `image_to_EIS` was **called with a path of image**, e. g. "images/geometry.png"
  and it means that z_file will have a form "images/geometry.z", i. e. changes only the extension to ".z"
- `save_also_image = ""
  - if "example.png" : the image will be saved with this name
  - if `!asZfile` : if `export_z_file != ""` and than the input image is copied with a name of `export_z_file` but with the extension ".png"
  - if `!input` : if `image_to_EIS` was **called with a path of image**
- `store_R_RC` : turns on the evaluation of R_RC element from two points of computed impedance and append the output to a specified file
  - if `= ""` : means *do nothing* (which is default)
  - if `= "storage.txt"` : append a line formated as <dateTtime> tab <input_file_name> tab <R_ohm> tab <R_pol> tab <C_pol> . If the function was called with a matrix (not an path to file), `"<matrix_input>"` is written instead of <input_file_name>.
  - if `= "storage.csv"` : csv extension works too
- `return_specific_impedance = true` : returns specific impedance in *ohms x cm*. Also values of fitted *R_RC* circuit will be returned in specific units.
- `L_el_mat = i_LSM` : defines material of the left electrode adjacent to the specified material matrix. Possible options are `i_LSM` (default) or `i_YSZ`.
- `R_el_mat = i_LSM` : the same for the right electrode material.

Advanced keyword parameters are 

- `iterative_solver = "auto"` : for small problems (under $15^3$ voxels) is used direct solver and iterative for larger ones.
  - if `= false` : the system of equations is solved by a direct LU solver (julian `\` operator)
  - if `= true` : the system is solved by iterative solver using Biconjugate gradient stabilized method with using Crout version of incomplete LU decomposition as a preconditioner.
- `fill_in_ratio = 12`: defines expected fill during the incomplete LU desomposition. It is a ratio of non-zero element count of aproximation of `L` matrix with respect to non-zero element count of original `A` matrix of the (sparse) linear system. Bigger `fill_in_ratio` means more realiable convergence of iterative solver but it requires more RAM space, which is the limiting factor for large systems (e.g. $10^6$ voxels and above).
- `tau = "auto"`: a drop criterion parameter in incomplete LU process. It is a absolute (not relative) treshold and elements of size under this treshold are forgotten. If `"auto"`, an adaptive process is run to estimate a feasible value for `tau` for which the defined `fill_in_ratio` is achieved.
- In a `par_study` dictionary, `"compute_tortuosity" => true` adds a `tor` column. This is a geometric shortest-path measure through the combined solid phase, not transport-derived tortuosity. Pore starting points are excluded and a domain with no solid path returns `Inf`.

#### Lower level API 
A user can directly access the computational core using function `material_matrix_to_impedance()` which accepts arguments these non-keyword arguments

- `material_matrix`
- `physical_parameters`

and keyword arguments

- `f_list`
- `iterative_solver`
- `verbose` : if `true`, print the frequency being solved
- `return_only_linsys` : if `true` : the linear system of equations is constructed and returned as a tuple `(A, b, current_measurement)` - matrix $A$ and right-hand-side $b$. The system is not solved.
- `tau`
- `fill_in_ratio`

### Real Example

```julialang
image_to_EIS(   [1 0 1; 0 1 2], 
                ["R_YSZ" => 73],
                #
                export_z_file="test.z", 
                return_R_RC=true,
                save_also_image="!asZfile",
                store_R_RC = ""
                )
```


## Generate random domain with structure

Automated generating of random structure is essential for statistical testing of system behavior. There are a few helping features using two main parameters. 

- `porosity` in [0, 1] : the ratio of pores over the total points (material points + pores) in the picture. 
- `LSM_ratio` in [0, 1] : probability that the material point will be LSM.

### Simple homogenous matrix
The simplest example is a homogenous domain of `dimensions = (m, n)`, where *m* is a number of rows and *n* a number of columns. Matrix of this type can be constructed via

```julialang
homogenous_matrix = generate_matrix(dimensions, porosity, LSM_ratio)
```

### Additional stuctural parameters
More specific structure can be defined using parameters which governs a tendency to group pixels together. Namely:

- `pore_cavitance = nothing`- a probability that next pore pixel will spawn next to already existing pore (if possible)
  - if `= nothing`, generating of random structure is fast because no "grouping algoritm" is involved.
  - if `= 0.0`, generating of random structure is slower because "grouping algoritm" is invoked
- `LSM_cavitance = 0.0`- a probability that next LSM pixel will spawn next to already existing LSM pixel (if possible). This keyword argument is working only if pore_cavitance keyword parameter is set to a real number (i.e. if `pore_cavitance = nothing`, then parameter `LSM_cavitance` has no effect to generated structure).

```julialang
pseudo_homogenous_matrix = generate_matrix(dimensions, porosity, LSM_ratio, 
                                            pore_cavitance = 0.4,
                                            LSM_cavitance = 0.2
                          )
```

### Rotation of domain
Generated domain can be rotated by function

```julialang
rotate_matrix!(domain, step, axis)
```

which rotates the input domain in place and returns it. Odd quarter-turns require equal lengths of the two rotated dimensions; otherwise use the nonmutating function, which returns an array with the corresponding dimensions exchanged. There is a version (withoug the exclamation mark `!`) which does not change the original domain and returns the new rotated domain

```julialang
new_domain = rotate_matrix(domain, step, axis)
```

The parameters are
- `step` specifing the number of simple rotations. 1 = 90°, 2 = 180°, 3 = 270°, 4 = 360°, 5 = 450°. And -1 = 90° backwards etc... 
- `axis` from the set `["x", "y", "z"]` specifies the axis of rotation. Note, the current flows along `"y"` axis, therefore rotation aroud this axis does not change the output impedance.

### Inspecting domain properties
Voxels in domain can be marked by a number, which will specify its material connectivity character. In particular,

 - 0 - pore
 - 1 - izolated material (not connected to electrodes)
 - 2 - material connected to left electrode
 - 3 - material connected to right electrode
 - 4 - material connected to both electrodes through the solid phase

 This investigation is dome by function 
 
 ```julialang
characterize_material_izolation(material_matrix)
 ```

 which returns matrix labeled as described above. In order to compute ratio of izolated material, there is a function

```julialang
izolated_material_ratio(material_matrix, wanted_ids = [1,2,3])
```
 
 which return the ratio of "wanted_material_points / all_material_points". The default value is `wanted_ids = [1]` for strongly izolated voxels.

### Structure using multiple submatrices
For more complicated domains composed of several different homogenous subdomains, there is a possibility to construct appropriate matrix. Suppose we want to construct *m* x *n* matrix consisting of 2 different submatrices. First, a list of submatrices must be created such that one submatrix is represented by its *location* (left upper corner and right lower corner) in the resulting matrix and *porosity* and *LSM_ratio*. 

```julialang
    # the structure for each submatrix in the list is >>
    # 
    # [left upper coord, right lower coord,  porosity,   LSM_ratio]
    
submatrix_list = [        
    [(1,1),             (10, 5),            0.2,          0.5],   # this is the first submatrix
    [(1,6),             (10, 10),           0.0,          1.0]    # this is the second submatrix
]

two_subdomain_matrix = generate_matrix(submatrix_list)
```

There can be more subdomains in `submatrix_list`. Dimensions of the resulting `two_subdomain_matrix` is computed as en rectangular envelope of all *locations* in `submatrix_list`. The submatrices can overlap (in this case, the latter has priority in evaluation of matrix), but **every pixel must be covered** by the submatrices.

#### Three column domain template
There is a template using the upper structure of defying submatrices, which generates a domain of 3 columns with defined material specification (*porosity* and *LSM_ratio*). In addition, there are contacts (of width 1 and optional height) on the left side providing an interesting distribution of electrical current through the system. The right side consists of a continuous one layer of LSM as a connection to conductive electrolyte. The obligate input parameters are 

- `LSM_ratio1`, `LSM_ratio2`, `LSM_ratio3`

Optional parameters (with default values) are 

- `porosity1=0.5, porosity2=0.5, porosity3=0.5`
- `positions_of_contacts=[15, 50]` : starting row for each LSM contact 
- `height_of_contacts=5`
- `column_width=5`
- `height=70`

in this default case, the LSM contacts will be between [15, 20] and [50, 55] pixels of the first column while the whole matrix has height of 70 pixels. The submatrix list can be then obtained by function `three_column_domain_template` and then inserted to `generate_matrix`.

```julialang
LSM_ratio1 = 0.2
LSM_ratio2 = 1.0
LSM_ratio3 = 0.5

template_submatrix_list = three_column_domain_template(LSM_ratio1, LSM_ratio2, LSM_ratio3,
                                     #                              
                                     column_width = 10,
                                     porosity1 = 0.0, porosity3 = 0.2,
                                     positions_of_contacts=[20, 45], height_of_contacts=4
                                 )
                        
three_column_matrix = generate_matrix(template_submatrix_list)
```

## Material matrix visualization

A material matrix can be saved to a file with specified `path` using `matrix_to_file` function. For example

```julialang
matrix_to_file("images/three_column_domain.png", three_column_matrix)
```

!["Three column domain"](images/three_column_domain.png?raw=true )

## Enlarge image

The images are often very small (50 pixels). So for a better view, there is a utility to make bigger resolution.

```julialang
enlarge_image(path; resolution=nothing, resize_factor=nothing, print_bool=false)
```

Parameters are

- `path`: path for the source image. The processed (enlarged) image is than saved to `"path_resized.XXX"`
- `resolution = nothing` you can choose if you want to specify the target resolution as an integer (it assumes squared pictures)
  - `resolution = 100` makes a resized picture with dimensions 100x100 pixels
- `resize_factor = nothing` or you can choose resizing factor, such that
  - `resize_factor = 3` means it will make 150x150 if the input image was 50x50
- `print_bool = false` ... if `true`, a name of resized picture is printed to terminal

## In-depth surface visualization

The function `make_shaded_view` takes a 3D structure `domain` and makes a 2D surface view such that it colors visible material pixels according to their depth in the structure. The basic usage is 

```julialang
make_shaded_view(domain, file_name, grad_depth=15, side=1, only_palette = false,
                            LSM_s = (0.3, 0.3, 0.3), LSM_e = (0.0, 0.0, 0.0),
                            YSZ_s = (1.0, 0.8, 0.0), YSZ_e = (0.5, 0.4, 0.0),
                            pall_func = x->x
                          )
```
Parameters are 

- `grad_depth = 15` 
- `only_palette = false`: if `true`, the output image is not the 2D surface view but only the palletes of colors
- `LSM_s`, `LSM_e`: starting and ending RGB color for gray LSM palette. Similarly, `YSZ_s` and `YSZ_e` for yellow YSZ pallete.
- `pall_func = x -> x`: user can define a function how pallete color should progress from start to end. The `x -> x` means linear progress. The input parameter `x` goes from 0 to 1 and the output of the function is expected also in the interval `[0, 1]`. There is a predefined symmetrical logaritmic function 
    - ` = sym_log(x, a)` which has the value `0.5` for `x = 0.5` and behaves as a logarithm in interval `[0.0, 0.5]` and is point-wise symmetric around the midpoint `x,y = (0.5, 0.5)`. The parameter `a` governs the "logarithmic steepness" in the beginng. `a = 1` is almost linear, `a = 1000` is very steep
- `side = 1`: defines, which side od the cube domain is inspected. In particular, the coding is 1 = do not rotate, 2 = 90° along x, 3 = 180° along x, 4 = 270° along x, 5 = 90° along y, 6 = -90° along y

## Matrix compilation of subimages

One can generate a lot of images using this package. There is also a tool to make a matrix of images based on row and column indexes:

```julialang
subimages_composition(folder_path = "data/set_of_images/", 
                      row_numbers = collect(0.2 : 0.1 : 0.7), 
                      col_numbers = collect(0.2 : 0.05 : 0.4), 
                      template = "img_por#(row)_LSM#(col).png",
)
```

where there parameters are

- `template = "img_por#(row)_LSM#(col)"` where the sequence `#(row)` will be replaced by a number from the list `row_numbers` nad the same holds for `#(col)`. Note that it is not `$(col)` so that it will not colide with the internal julia string interpolation.
- `background_color = (1,1,1)` RGB definition of background color. White is default.
- `del_width = 10` number of pixels separating subimages
- `out_filename = "compilation.png"` and output file is located in `folder_path`


### Connectivity and bounded generation

`generate_matrix` supports the following keywords for both random domains and lists of subdomains:

- `check_connectivity=true`: retain only packings with a solid path between the two electrodes. This samples a distribution **conditioned on connectivity**; it is not an unconditional percolation study.
- `check_connectivity=false`: return the first generated packing, including disconnected realizations. Use this for unbiased sampling of the generator and report the connected fraction separately.
- `max_attempts=1000`: maximum total number of generated candidates. Failure throws an error, rather than returning a scalar in place of an array. The legacy `recursion_depth` keyword remains accepted as an alias for this total attempt limit.
- `generation_stats=nothing`: optionally pass a `Dict{String,Any}()` to receive `"attempts"`, `"rejected"`, and `"connected"`. On failure the dictionary contains the final attempt count. These are statistics of this call; aggregate them across a study as needed.

For composite structures each subdomain is generated once without a connectivity restriction. Contacts and overlaps are applied before testing the **assembled** packing. Retrying regenerates the whole candidate. This permits intentionally pore-only subdomains.

```julia
stats = Dict{String,Any}()
packing = generate_matrix((30, 30), 0.5, 0.5;
    check_connectivity=false, generation_stats=stats)
classification = characterize_material_izolation(packing)
isolated = classification .== 1
```

Connectivity, isolation classification, and shortest paths currently treat LSM and YSZ together as one solid phase. They do not separately assess ionic and electronic percolation. Pores in the electrical network retain the finite `R_pore`, so solid connectivity is a geometric classification, not a guarantee of exactly zero electrical leakage.

### Numerical validation

Both direct and iterative results are checked for finite potentials, the true residual `norm(A*x-b)`, and agreement of left and right electrode currents. Invalid results raise an error instead of silently becoming an impedance. Both `material_matrix_to_impedance` and `image_to_EIS` accept:

- `solver_rtol=1e-8`, `solver_atol=0.0`: true residual must be at most `solver_atol + solver_rtol*norm(b)`.
- `current_rtol=1e-6`, `current_atol=1e-12`: current imbalance must be at most `current_atol + current_rtol*max(abs(I_left), abs(I_right))`. Absolute current tolerance is in raw network current units, before specific-impedance normalization.
- `max_mv_products=2000`: BiCGStab matrix-vector product limit. If BiCGStab breaks down or fails the residual check, a restarted GMRES retry has at most this many iterations; it does not fall back to a direct factorization.
- `solver_diagnostics=nothing`: optionally pass an empty vector to collect frequency, solver method, residual, both electrode currents and their absolute difference. Entries correspond to directly solved frequencies, including only the two input frequencies in TPE mode.

The automatic ILU tolerance search uses at most 40 trial factorizations of its sample matrix. An unattainable fill target does not cause unlimited iteration. The fill target is an estimate on that sample, not a strict memory bound for the full matrix.

### Running the tests

The GitHub Actions workflow tests Julia 1.8 and 1.10 using Python 3.10, the recorded Julia dependencies and Matplotlib 3.8.4 with NumPy 1.26.4. The Matplotlib pin is needed by the recorded PyPlot version (newer Matplotlib releases remove an API it uses). Setuptools 75.8.0 supplies the legacy `distutils` module when using Python 3.12 locally. Install `test/python-requirements.txt` in an isolated Python environment before running:

```sh
MPLBACKEND=Agg PYTHON=/path/to/python julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
MPLBACKEND=Agg julia --project -e 'using Test; include("test/runtests.jl")'
```

Tests include analytical 2D/3D spectra, LSM and YSZ contacts, exact frequency endpoints, solver agreement, invalid-solution rejection, bounded generation, rotations, isolation, and geometric shortest paths.

## Acknowledgement

This work was supported by the German Research Foundation, DFG project no. FU 316/14-1, and by the Czech Science Foundation, GAČR project no. 19-14244J.












