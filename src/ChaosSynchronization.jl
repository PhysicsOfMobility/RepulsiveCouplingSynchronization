module ChaosSynchronization

include("ChaosSynchronizationODEs.jl")
export 
    SystemParameters,
    roesslerODE, 
    roesslerODE_coupled, 
    roesslerODE_coupled_linear,
    roesslerODE_coupled_linear_focal_test,
    lorenzODE,
    lorenzODE_coupled,
    lorenzODE_coupled_linear,
    chenODE,
    chenODE_coupled,
    chenODE_coupled_linear,
    chuaODE,
    chuaODE_coupled,
    chuaODE_coupled_linear,
    roessler_backward_linear,
    SystemParameters_discontinuity_detection,
    roesslerODE_coupled_discontinuity_detection,
    roesslerODE_coupled_linear_discontinuity_detection,
    coupling_condition,
    #coupling_affect!,
    #coupling_init,
    coupling_callback,
    nonlinear_coupling_condition,
    #nonlinear_coupling_init,
    nonlinear_coupling_callback
    


include("ChaosSynchronizationHelper.jl")
export  
    setup_system, 
    generate_attractor_sample, 
    generate_attractor_plot, 
    generate_attractortiling_uniform, 
    generate_attractortiling_equidistant, 
    generate_couplingregions, 
    generate_couplingregions_reshuffled, 
    simulate_nonlinear_trajectory,
    simulate_linearized_trajectory,
    compute_maximum_lyapunov_exponent_linearized,
    scan_maximum_lyapunov_exponent,
    compute_coupling_region_lyap_effect,
    find_next_optimal_coupling_region,
    simulate_nonlinear_optimal_coupling_trajectory,
    simulate_linearized_optimal_coupling_trajectory,
    simple_test__compute_maximum_lyapunov_exponent_linearized,
    simple_test__scan_maximum_lyapunov_exponent,
    compute_coupling_fraction,
    is_in_coupling_region,
    compute_stable_direction_uncoupled_linearized,
    simulate_synch_time_distribution,
    compute_maximum_lyapunov_exponent_linearized_discontinuity_handling,
    find_next_optimal_coupling_region_discontinuity_handling,
    simulate_nonlinear_trajectory_discontinuity_handling,
    simulate_linearized_trajectory_discontinuity_handling,
    simulate_nonlinear_optimal_coupling_trajectory_discontinuity_handling,
    simulate_linearized_optimal_coupling_trajectory_discontinuity_handling,
    scan_maximum_lyapunov_exponent_discontinuity_handling


include("ChaosSynchronizationMakieStyle.jl")

include("ChaosSynchronizationPlot.jl")
export 
    attractor_plot,
    couplingregion_plot,
    couplingregion_plot_all,
    optimal_coupling_regions_plot,
    attractor_difference_direction_plot,
    attractor_difference_direction_with_optimal_coupling_regions_plot,
    couplingregion_lyap_effect_plot,
    synch_error_plot_nonlinear,
    synch_error_plot_linearized,
    lyap_scan_plot

include("ChaosSynchronizationOptimization.jl")
export 
    optimize_coupling_regions,
    optimize_coupling_regions_discontinuity_handling

end