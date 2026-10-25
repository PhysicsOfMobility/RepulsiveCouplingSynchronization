using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays
import Printf
import DelimitedFiles

coupling_strength = StaticArrays.@SVector[-1.0, 0.0, 0.0]

number_of_coupling_regions_to_plot = 15

coupling_regions_filename = "couplingregions_equidist__dist_0.50__cfrac_0.01.tsv"


optimal_coupling_regions_filename = string(coupling_regions_filename[1:end-4], "__optimal_coupling__alpha_",Printf.@sprintf("%.2f",coupling_strength[1]),"_",Printf.@sprintf("%.2f",coupling_strength[2]),"_",Printf.@sprintf("%.2f",coupling_strength[3]),".tsv")

simulate_synch_time_distribution(roesslerODE_coupled, "roessler", optimal_coupling_regions_filename, coupling_strength)
