using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays 

#generate_attractortiling_uniform("roessler", attractor_fraction_tiling_list = [0.001, 0.002])
#generate_attractortiling_equidistant("roessler", min_dist_tiling_list = [0.5], shuffle_attractor_points = true)

generate_couplingregions_reshuffled("roessler", attractor_fraction_coupling_list = [0.01, 0.02, 0.05], min_dist_tiling_list = [0.5], attractor_fraction_tiling_list = Float64[])


