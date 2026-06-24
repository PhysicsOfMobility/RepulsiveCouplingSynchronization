using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays 

setup_system(
    lorenzODE, "lorenz"; 
    x0 = StaticArrays.SVector{3,Float64}(1.0,1.0,1.0),
    attractor_fraction_coupling_list = Float64[0.02], 
    min_dist_tiling_list = Float64[1.0], 
    attractor_fraction_tiling_list = Float64[]
    )

