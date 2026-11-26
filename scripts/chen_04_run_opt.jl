using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays
import Printf


coupling_strength = StaticArrays.@SVector[1.0, 0.0, 0.0]

number_of_coupling_regions_to_compute = 20

coupling_regions_filename = "couplingregions_equidist__dist_1.50__cfrac_0.01.tsv"

optimize_coupling_regions(
    chenODE_coupled_linear, 
    "chen", 
    coupling_regions_filename, 
    coupling_strength;
    number_of_coupling_regions = number_of_coupling_regions_to_compute, 
    show_progress = true
)

optimal_coupling_region_filename = string(coupling_regions_filename[1:end-4], "__optimal_coupling__alpha_",Printf.@sprintf("%.2f",coupling_strength[1]),"_",Printf.@sprintf("%.2f",coupling_strength[2]),"_",Printf.@sprintf("%.2f",coupling_strength[3]),".tsv")

compute_coupling_fraction(
    "chen", 
    optimal_coupling_region_filename
)