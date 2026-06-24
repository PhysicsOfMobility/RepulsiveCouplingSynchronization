using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays 

setup_system(
    roesslerODE, "roessler"; 
    attractor_fraction_coupling_list = Float64[0.01,0.05], 
    min_dist_tiling_list = Float64[0.5], 
    attractor_fraction_tiling_list = Float64[]
    )

