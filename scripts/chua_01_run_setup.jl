using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays 

setup_system(
    chuaODE, "chua"; 
    x0 = StaticArrays.SVector{3,Float64}(1.0,1.0,1.0),
    attractor_fraction_coupling_list = Float64[0.01], 
    min_dist_tiling_list = Float64[0.2], 
    attractor_fraction_tiling_list = Float64[]
    )

