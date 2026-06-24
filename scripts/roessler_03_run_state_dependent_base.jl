using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays
import Printf

coupling_strength = StaticArrays.@SVector[1.0, 0.0, 0.0]

coupling_region_filename_list = 
[
    "couplingregions_equidist__dist_0.50__cfrac_0.01.tsv", 
]

for coupling_region_filename in coupling_region_filename_list
    compute_coupling_region_lyap_effect(
        roesslerODE_coupled_linear, 
        "roessler", 
        coupling_region_filename, 
        coupling_strength,
        show_progress = true
    )
end

#for coupling_region_filename in coupling_region_filename_list
#    data_filename = string(coupling_region_filename[1:end-4], "__lyap_effect__alpha_",Printf.@sprintf("%.2f",coupling_strength[1]),"_",Printf.@sprintf("%.2f",coupling_strength[2]),"_",Printf.@sprintf("%.2f",coupling_strength[3]),".tsv")
#    couplingregion_lyap_effect_plot(
#        "roessler", 
#        data_filename
#    )
#end

