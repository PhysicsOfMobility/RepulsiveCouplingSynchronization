import StaticArrays 
import DifferentialEquations 
import DelimitedFiles 
import DynamicalSystems
import Distances
import Statistics
import Random
import Printf
import DrWatson

using ProgressMeter

function optimize_coupling_regions(
    ode::Function, 
    system_name::AbstractString, 
    coupling_region_filename::String, 
    coupling_strength::StaticArrays.SVector{3, Float64}; 
    number_of_coupling_regions::Int64 = 10, 
    number_of_coupling_region_samples::Int64 = -1, 
    initial_state_index::Int64 = 1, 
    fixed_x0::StaticArrays.SVector{7,Float64} = StaticArrays.SVector{7,Float64}(0.0,0.0,0.0,0.5,0.5,0.714,0.0), 
    delta_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(1.0, 1.0, 1.0), 
    simulation_time::Float64 = 100000.0, 
    transient_time::Float64 = 1000.0, 
    number_of_lyap_samples::Int64 = 10,
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4,
    show_progress::Bool = false
    )

    best_coupling_regions = Matrix{Float64}(undef, number_of_coupling_regions, 8)

    #start without any coupling
    systemparameters_coupled = SystemParameters(coupling_strength, 0, StaticArrays.SVector{3, Float64}[], Float64[])

    for k = 1 : number_of_coupling_regions
        print("finding region ",k," out of ",number_of_coupling_regions,"\n")
        best_next_coupling_region = find_next_optimal_coupling_region(
            ode, 
            system_name, 
            systemparameters_coupled,
            coupling_region_filename; 
            number_of_coupling_region_samples = number_of_coupling_region_samples,
            initial_state_index = initial_state_index,
            fixed_x0 = fixed_x0,
            delta_x0 = delta_x0,
            simulation_time = simulation_time,
            transient_time = transient_time,
            number_of_lyap_samples = number_of_lyap_samples,
            max_integration_timestep = max_integration_timestep,
            absolute_tolerance = absolute_tolerance,
            relative_tolerance = relative_tolerance,
            show_progress = show_progress
        )

        best_coupling_regions[k,:] = best_next_coupling_region

        best_next_coupling_region_center = StaticArrays.SVector{3, Float64}(best_next_coupling_region[1:3])
        best_next_coupling_region_radius = best_next_coupling_region[4]
        systemparameters_coupled = SystemParameters(systemparameters_coupled.couplingstrength, systemparameters_coupled.numberofcouplingregions + 1, vcat(systemparameters_coupled.couplingpoints,[best_next_coupling_region_center]), vcat(systemparameters_coupled.couplingradii,[best_next_coupling_region_radius]))
    end

    output_filename = string(coupling_region_filename[1:end-4], "__optimal_coupling__alpha_",Printf.@sprintf("%.2f",coupling_strength[1]),"_",Printf.@sprintf("%.2f",coupling_strength[2]),"_",Printf.@sprintf("%.2f",coupling_strength[3]),".tsv")
    DelimitedFiles.writedlm(DrWatson.datadir(system_name, output_filename), best_coupling_regions, '\t')

end
