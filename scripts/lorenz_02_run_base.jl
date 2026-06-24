using DrWatson
@quickactivate :ChaosSynchronization

import StaticArrays 

#systemparameters_uncoupled = SystemParameters(
#    StaticArrays.@SVector[0.0, 0.0, 0.0],
#    0, 
#    StaticArrays.SVector{3, Float64}[], 
#    Float64[]
#)

#simulate_nonlinear_trajectory(lorenzODE_coupled,"lorenz",systemparameters_uncoupled; sampling = true, output_filename = "nonlinear__uncoupled_trajectory_plot.tsv", sampling_timestep = 0.01)
#simulate_linearized_trajectory(lorenzODE_coupled_linear,"lorenz",systemparameters_uncoupled; sampling = true, output_filename = "linearized__uncoupled_trajectory_plot.tsv", sampling_timestep = 0.01)
#simulate_linearized_trajectory(lorenzODE_coupled_linear,"lorenz",systemparameters_uncoupled; sampling = true, output_filename = "linearized__uncoupled_trajectory_sample.tsv", sampling_timestep = 0.1)

#attractor_difference_direction_plot("lorenz", difference_filename_ovewrite = "linearized__uncoupled_trajectory_sample.tsv")

#synch_error_plot_nonlinear("lorenz", "nonlinear__uncoupled_trajectory_plot.tsv")
#synch_error_plot_linearized("lorenz", "linearized__uncoupled_trajectory_plot.tsv")

systemparameters_fullycoupled = SystemParameters(
    StaticArrays.@SVector[1.0, 0.0, 0.0], 
    1, 
    [StaticArrays.@SVector[0.0, 0.0, 0.0]], 
    [1000.0]
)

#simulate_nonlinear_trajectory(lorenzODE_coupled,"lorenz",systemparameters_fullycoupled; sampling = true, output_filename = "nonlinear__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00.tsv", sampling_timestep = 0.01)
#simulate_linearized_trajectory(lorenzODE_coupled_linear,"lorenz",systemparameters_fullycoupled; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00.tsv", sampling_timestep = 0.01)
#simulate_linearized_trajectory(lorenzODE_coupled_linear,"lorenz",systemparameters_fullycoupled; sampling = true, output_filename = "linearized__xcoupled_trajectory_sample__alpha_1.00_0.00_0.00.tsv", sampling_timestep = 0.1)

#attractor_difference_direction_plot("lorenz", difference_filename_ovewrite = "linearized__xcoupled_trajectory_sample__alpha_1.00_0.00_0.00.tsv")

#synch_error_plot_nonlinear("lorenz", "nonlinear__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00.tsv")
#synch_error_plot_linearized("lorenz", "linearized__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00.tsv")

coupling_strength_range = range(0.0, 10.0, step=0.01)

scan_maximum_lyapunov_exponent(
    lorenzODE_coupled_linear, 
    "lorenz", 
    systemparameters_fullycoupled; 
    coupling_strength_range = [coupling_strength_range, zeros(size(coupling_strength_range,1)), zeros(size(coupling_strength_range,1))],
    output_filename_ovewrite = "lyap_scan__fullycoupled__alpha_1.00_0.00_0.00.tsv",
    show_progress = true
    )

#lyap_scan_plot("lorenz", "lyap_scan__fullycoupled__alpha_1.00_0.00_0.00.tsv")