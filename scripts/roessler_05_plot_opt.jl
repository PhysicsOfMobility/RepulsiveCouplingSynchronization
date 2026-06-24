using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays
import Printf
import DelimitedFiles

coupling_strength = StaticArrays.@SVector[-1.0, 0.0, 0.0]

number_of_coupling_regions_to_plot = 15

coupling_regions_filename = "couplingregions_equidist__dist_0.50__cfrac_0.01.tsv"


optimal_coupling_regions_filename = string(coupling_regions_filename[1:end-4], "__optimal_coupling__alpha_",Printf.@sprintf("%.2f",coupling_strength[1]),"_",Printf.@sprintf("%.2f",coupling_strength[2]),"_",Printf.@sprintf("%.2f",coupling_strength[3]),".tsv")
#optimal_coupling_regions_plot("roessler",optimal_coupling_regions_filename; number_of_coupling_regions = number_of_coupling_regions_to_plot)

simulate_nonlinear_optimal_coupling_trajectory(
    roesslerODE_coupled,
    "roessler",
    optimal_coupling_regions_filename,
    coupling_strength;
    #coupling_strength_sim;
    number_of_coupling_regions = number_of_coupling_regions_to_plot,
    sampling = true,
    sampling_timestep = 10.0,
    transient_time = 4500.0,
    simulation_time = 1.0e6,
    output_suffix = "plot"
)



systemparameters_uncoupled = SystemParameters(
    StaticArrays.@SVector[0.0, 0.0, 0.0],
    0, 
    StaticArrays.SVector{3, Float64}[], 
    Float64[]
)

simulate_nonlinear_trajectory(
    roesslerODE_coupled,
    "roessler",
    systemparameters_uncoupled; 
    sampling = true, 
    output_filename = "nonlinear__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", 
    sampling_timestep = 10.0, 
    transient_time = 4500.0, 
    simulation_time = 1e6
)



simulate_linearized_optimal_coupling_trajectory(
    roesslerODE_coupled_linear,
    "roessler",
    optimal_coupling_regions_filename,
    coupling_strength;
    #coupling_strength_sim;
    number_of_coupling_regions = number_of_coupling_regions_to_plot,
    sampling = true,
    transient_time = 4560.0,
    simulation_time = 10000.0,
    sampling_timestep = 0.1,
    output_suffix = "plot",
)

#simulate_linearized_optimal_coupling_trajectory(
#    roesslerODE_coupled_linear,
#    "roessler",
#    optimal_coupling_regions_filename,
#    coupling_strength;
#    #coupling_strength_sim;
#    number_of_coupling_regions = number_of_coupling_regions_to_plot,
#    sampling = true,
#    transient_time = 4560.0,
#    sampling_timestep = 0.1,
#    output_suffix = "sample"
#)

#attractor_difference_plot_filename = string("linearized__",optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_sample.tsv")
#attractor_difference_plot_output_filename = string("linearized__",optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_sample.png")
#attractor_difference_direction_plot("roessler", difference_filename_ovewrite = attractor_difference_plot_filename, output_filename_overwrite = attractor_difference_plot_output_filename)
#attractor_difference_direction_with_optimal_coupling_regions_plot(
#    "roessler", 
#    optimal_coupling_regions_filename, 
#    attractor_difference_plot_filename; 
#    number_of_coupling_regions = number_of_coupling_regions_to_plot
#)
#
#synch_error_plot_nonlinear_filename = string("nonlinear__", optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_plot.tsv")
#synch_error_plot_nonlinear("roessler", synch_error_plot_nonlinear_filename)
#synch_error_plot_nonlinear_filename = string("linearized__", optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_plot.tsv")
#synch_error_plot_linearized("roessler", synch_error_plot_nonlinear_filename)



#optimal_coupling_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler", optimal_coupling_regions_filename), '\t', Float64)
#
#number_of_optimal_coupling_points = number_of_coupling_regions_to_plot
#optimal_coupling_points = [StaticArrays.SVector{3, Float64}(optimal_coupling_data[i, 1:3]) for i in 1:number_of_optimal_coupling_points]
#optimal_coupling_radii = optimal_coupling_data[:, 4]
#
#optimal_coupling_parameters = SystemParameters(coupling_strength, number_of_optimal_coupling_points, optimal_coupling_points, optimal_coupling_radii)
#
#coupling_strength_range = range(-10.0, 10.0, step=0.01)
#
#scan_maximum_lyapunov_exponent(
#    roesslerODE_coupled_linear, 
#    "roessler", 
#    optimal_coupling_parameters; 
#    coupling_strength_range = [coupling_strength_range, zeros(size(coupling_strength_range,1)), zeros(size(coupling_strength_range,1))],
#    output_filename_ovewrite = "lyap_scan__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_1.00_0.00_0.00_ncr_15.tsv",
#    show_progress = true
#    )
#
#lyap_scan_plot("roessler", "lyap_scan__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_1.00_0.00_0.00_ncr_15.tsv")