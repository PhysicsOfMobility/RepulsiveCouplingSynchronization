using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays
import Printf
import DelimitedFiles

coupling_strength = StaticArrays.@SVector[5.0, 0.0, 0.0]
#coupling_strength_sim = StaticArrays.@SVector[10.0, 0.0, 0.0]

number_of_coupling_regions_to_plot = 10

coupling_regions_filename = "couplingregions_equidist__dist_0.50__cfrac_0.05.tsv"


optimal_coupling_regions_filename = string(coupling_regions_filename[1:end-4], "__optimal_coupling__alpha_",Printf.@sprintf("%.2f",coupling_strength[1]),"_",Printf.@sprintf("%.2f",coupling_strength[2]),"_",Printf.@sprintf("%.2f",coupling_strength[3]),".tsv")

optimal_coupling_regions_plot("lorenz",optimal_coupling_regions_filename; number_of_coupling_regions = number_of_coupling_regions_to_plot)

simulate_nonlinear_optimal_coupling_trajectory(
    lorenzODE_coupled,
    "lorenz",
    optimal_coupling_regions_filename,
    coupling_strength;
    #coupling_strength_sim;
    number_of_coupling_regions = number_of_coupling_regions_to_plot,
    sampling = true,
    #sampling_timestep = 1.0,
    sampling_timestep = 0.05,
    #transient_time = 4500.0,
    #simulation_time = 1.0e6,
    simulation_time = 1000.0,
    output_suffix = "plot"
)

simulate_linearized_optimal_coupling_trajectory(
    lorenzODE_coupled_linear,
    "lorenz",
    optimal_coupling_regions_filename,
    coupling_strength;
    #coupling_strength_sim;
    number_of_coupling_regions = number_of_coupling_regions_to_plot,
    sampling = true,
    #transient_time = 4550.0,
    simulation_time = 10000.0,
    sampling_timestep = 0.01,
    output_suffix = "plot",
)

simulate_linearized_optimal_coupling_trajectory(
    lorenzODE_coupled_linear,
    "lorenz",
    optimal_coupling_regions_filename,
    coupling_strength;
    #coupling_strength_sim;
    number_of_coupling_regions = number_of_coupling_regions_to_plot,
    sampling = true,
    #transient_time = 4550.0,
    sampling_timestep = 0.1,
    output_suffix = "sample"
)

attractor_difference_plot_filename = string("linearized__",optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_sample.tsv")
attractor_difference_plot_output_filename = string("linearized__",optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_sample.png")
attractor_difference_direction_plot("lorenz", difference_filename_ovewrite = attractor_difference_plot_filename, output_filename_overwrite = attractor_difference_plot_output_filename)
attractor_difference_direction_with_optimal_coupling_regions_plot(
    "lorenz", 
    optimal_coupling_regions_filename, 
    attractor_difference_plot_filename; 
    number_of_coupling_regions = number_of_coupling_regions_to_plot
)

synch_error_plot_nonlinear_filename = string("nonlinear__", optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_plot.tsv")
synch_error_plot_nonlinear("lorenz", synch_error_plot_nonlinear_filename)
synch_error_plot_nonlinear_filename = string("linearized__", optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_coupling_regions_to_plot,"__trajectory_plot.tsv")
synch_error_plot_linearized("lorenz", synch_error_plot_nonlinear_filename)



#
#optimal_coupling_data = DelimitedFiles.readdlm(DrWatson.datadir("lorenz", optimal_coupling_regions_filename), '\t', Float64)
#
#number_of_optimal_coupling_points = number_of_coupling_regions_to_plot
#optimal_coupling_points = [StaticArrays.SVector{3, Float64}(optimal_coupling_data[i, 1:3]) for i in 1:number_of_optimal_coupling_points]
#optimal_coupling_radii = optimal_coupling_data[:, 4]
#
#optimal_coupling_parameters = SystemParameters(coupling_strength, number_of_optimal_coupling_points, optimal_coupling_points, optimal_coupling_radii)
#
#coupling_strength_range = range(0.0, 10.0, step=0.01)
#
#scan_maximum_lyapunov_exponent(
#    lorenzODE_coupled_linear, 
#    "lorenz", 
#    optimal_coupling_parameters; 
#    coupling_strength_range = [coupling_strength_range, zeros(size(coupling_strength_range,1)), zeros(size(coupling_strength_range,1))],
#    output_filename_ovewrite = "lyap_scan__couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00_ncr_10.tsv",
#    show_progress = true
#    )
#
#lyap_scan_plot("lorenz", "lyap_scan__couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00_ncr_10.tsv")