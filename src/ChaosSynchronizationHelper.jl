import StaticArrays 
import DifferentialEquations 
import DelimitedFiles 
import DynamicalSystems
import SimpleDiffEq
import Distances
import Statistics
import Random
import Printf
import DrWatson

using ProgressMeter

function generate_attractor_sample(
    ode::Function, 
    system_name::AbstractString; 
    x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.0,0.0,0.0), 
    sampling_timestep_sample::Float64 = 1.0, 
    simulation_time_sample::Float64 = 100000.0, 
    transient_time::Float64 = 1000.0, 
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4 
    )

    export_filename::String = "attractor_sample.tsv"

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            ode, 
            x0,
            diffeq = (
                save_everystep = false, 
                save_start = false,
                dtmax = max_integration_timestep, 
                abstol = absolute_tolerance,
                reltol = relative_tolerance
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)
    x0 = DynamicalSystems.current_state(ds)
    DynamicalSystems.reinit!(ds, x0, t0 = 0.0)

    u, t = DynamicalSystems.trajectory(ds, simulation_time_sample; Δt = sampling_timestep_sample)
    u = u[Random.randperm(size(u, 1)),:]

    DelimitedFiles.writedlm(DrWatson.datadir(system_name, export_filename), u, '\t')

    nothing
end

function generate_attractor_plot(
    ode::Function, 
    system_name::AbstractString; 
    x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.0,0.0,0.0), 
    sampling_timestep::Float64 = 0.01, 
    simulation_time::Float64 = 2000.0, 
    transient_time::Float64 = 1000.0, 
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4 
    )
    export_filename::String = "attractor_plot.tsv"

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            ode, 
            x0,
            diffeq = (
                save_everystep = false, 
                save_start = false, 
                dtmax = max_integration_timestep,
                abstol = absolute_tolerance,
                reltol = relative_tolerance
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)
    x0 = DynamicalSystems.current_state(ds)
    DynamicalSystems.reinit!(ds, x0, t0 = 0.0)

    u, t = DynamicalSystems.trajectory(ds, simulation_time; Δt = sampling_timestep)

    DelimitedFiles.writedlm(DrWatson.datadir(system_name,export_filename), hcat(t,u), '\t')

    nothing
end

function generate_attractortiling_equidistant(
    system_name::AbstractString; 
    min_dist_tiling_list::Vector{Float64} = [0.25, 0.5]
    )

    attractor_filename::String = "attractor_sample.tsv"
    attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
    N = size(attractorpoints, 1)
    
    for min_dist in min_dist_tiling_list

        min_dist_sq = min_dist*min_dist
        tiling_filename::String = string("attractortiling_equidist__dist_",Printf.@sprintf("%.2f",min_dist),".tsv")

        selected_indices = Int[]
        remaining_indices = Set([i for i in 1:N])

        current_index = 1
        push!(selected_indices, current_index)
        delete!(remaining_indices, current_index)

        min_dist_to_selected_sq = [Inf for i in 1:N]

        npoints = 1

        while !isempty(remaining_indices)
            for i in remaining_indices
                min_dist_to_selected_sq[i] = min(min_dist_to_selected_sq[i], Distances.SqEuclidean()(attractorpoints[i,:],attractorpoints[current_index,:]))
                if(min_dist_to_selected_sq[i] < min_dist_sq)
                    delete!(remaining_indices, i)
                end
            end

            if(isempty(remaining_indices))
                break
            end

            candidate_idx = -1
            min_min_dist_sq = Inf

            for idx in remaining_indices
                if min_dist_to_selected_sq[idx] < min_min_dist_sq
                    min_min_dist_sq = min_dist_to_selected_sq[idx]
                    candidate_idx = idx
                end
            end

            current_index = candidate_idx
            push!(selected_indices, current_index)
            delete!(remaining_indices, current_index)
            npoints += 1
        end

        couplingpoints = []
        for idx in selected_indices
            push!(couplingpoints, attractorpoints[idx, :])
        end

        DelimitedFiles.writedlm(DrWatson.datadir(system_name,tiling_filename), couplingpoints, '\t')
    end

    nothing
end

function generate_attractortiling_uniform(
    system_name::AbstractString; 
    attractor_fraction_tiling_list::Vector{Float64} = [0.001, 0.002]
    )

    attractor_filename::String = "attractor_sample.tsv"
    attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
    N = size(attractorpoints, 1)
    
    for attractor_fraction_tiling in attractor_fraction_tiling_list
        attractor_fraction_index = ceil(Int, N * attractor_fraction_tiling)
        attractor_fraction_index = clamp(attractor_fraction_index, 1, N)

        tiling_filename::String = string("attractortiling_uniform__frac_",Printf.@sprintf("%.3f",attractor_fraction_tiling),".tsv")

        selected_indices = Int64[]
        remaining_indices = Set([i for i in 1:N])
        current_attractorfraction_sqdistance = 0.0
        current_attractorfraction_distance = 0.0

        current_sqdistance_list = zeros(N)
        current_distance_list = zeros(N)

        current_index = 1
        push!(selected_indices, current_index)
        delete!(remaining_indices, current_index)

        min_dist_to_selected = [Inf for i in 1:N]

        npoints = 1

        while !isempty(remaining_indices)

            for i in 1:N
                current_sqdistance_list[i] = Distances.SqEuclidean()(attractorpoints[i,:],attractorpoints[current_index,:])
                current_distance_list[i] = sqrt(current_sqdistance_list[i])
            end

            current_attractorfraction_sqdistance = partialsort(current_sqdistance_list, attractor_fraction_index)
            current_attractorfraction_distance = sqrt(current_attractorfraction_sqdistance)

            for i in remaining_indices
                if(current_sqdistance_list[i] < current_attractorfraction_sqdistance)
                    delete!(remaining_indices, i)
                else
                    min_dist_to_selected[i] = min(min_dist_to_selected[i], current_distance_list[i] - current_attractorfraction_distance)
                end
            end

            if(isempty(remaining_indices))
                break
            end

            candidate_idx = -1
            min_min_dist = Inf

            for idx in remaining_indices
                if min_dist_to_selected[idx] < min_min_dist
                    min_min_dist = min_dist_to_selected[idx]
                    candidate_idx = idx
                end
            end

            current_index = candidate_idx
            push!(selected_indices, current_index)
            delete!(remaining_indices, current_index)
            npoints += 1
        end

        couplingpoints = []
        for idx in selected_indices
            push!(couplingpoints, attractorpoints[idx, :])
        end

        DelimitedFiles.writedlm(DrWatson.datadir(system_name,tiling_filename), couplingpoints, '\t')
    end

    nothing
end

function generate_couplingregions(
    system_name::AbstractString; 
    attractor_fraction_coupling_list::Vector{Float64} = [0.01, 0.02, 0.05], 
    min_dist_tiling_list::Vector{Float64} = [0.25, 0.5], 
    attractor_fraction_tiling_list::Vector{Float64} = [0.001, 0.002]
    )

    attractor_filename::String = string("attractor_sample.tsv")
    attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
    N = size(attractorpoints, 1)
    
    for min_dist_tiling in min_dist_tiling_list
        tiling_filename::String = string("attractortiling_equidist__dist_",Printf.@sprintf("%.2f",min_dist_tiling),".tsv")
        couplingpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,tiling_filename), '\t', Float64)

        ncoupling = size(couplingpoints, 1)
        nattractor = size(attractorpoints, 1)
        couplingRadius = zeros(ncoupling,size(attractor_fraction_coupling_list,1))
        couplingIndices = [ clamp(ceil(Int, nattractor * attractor_fraction_coupling_list[k]),1,nattractor) for k in 1:size(attractor_fraction_coupling_list,1)]

        Threads.@threads for i in 1:ncoupling
            distances_sq = Vector{Float64}(undef, nattractor)

            for j in 1:nattractor
                distances_sq[j] = Distances.SqEuclidean()(couplingpoints[i, :], attractorpoints[j, :])
            end

            for k in 1:size(attractor_fraction_coupling_list,1)
                couplingRadius[i,k] = partialsort(distances_sq, couplingIndices[k])
            end
        end

        for k in 1:size(attractor_fraction_coupling_list,1)
            couplingregion_filename::String = string("couplingregions_equidist__dist_",Printf.@sprintf("%.2f",min_dist_tiling),"__cfrac_",Printf.@sprintf("%.2f",attractor_fraction_coupling_list[k]),".tsv")
            DelimitedFiles.writedlm(DrWatson.datadir(system_name,couplingregion_filename), hcat(couplingpoints,sqrt.(couplingRadius[:,k]),couplingRadius[:,k]), '\t')
        end
    end


    for attractor_fraction_tiling in attractor_fraction_tiling_list
        tiling_filename = string("attractortiling_uniform__frac_",Printf.@sprintf("%.3f",attractor_fraction_tiling),".tsv")
        couplingpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,tiling_filename), '\t', Float64)

        ncoupling = size(couplingpoints, 1)
        nattractor = size(attractorpoints, 1)
        couplingRadius = zeros(ncoupling,size(attractor_fraction_coupling_list,1))
        couplingIndices = [ clamp(ceil(Int, nattractor * attractor_fraction_coupling_list[k]),1,nattractor) for k in 1:size(attractor_fraction_coupling_list,1)]

        Threads.@threads for i in 1:ncoupling
            distances_sq = Vector{Float64}(undef, nattractor)

            for j in 1:nattractor
                distances_sq[j] = Distances.SqEuclidean()(couplingpoints[i, :], attractorpoints[j, :])
            end

            for k in 1:size(attractor_fraction_coupling_list,1)
                couplingRadius[i,k] = partialsort(distances_sq, couplingIndices[k])
            end
        end

        for k in 1:size(attractor_fraction_coupling_list,1)
            couplingregion_filename::String = string("couplingregions_uniform__frac_",Printf.@sprintf("%.3f",attractor_fraction_tiling),"__cfrac_",Printf.@sprintf("%.2f",attractor_fraction_coupling_list[k]),".tsv")
            DelimitedFiles.writedlm(DrWatson.datadir(system_name, couplingregion_filename), hcat(couplingpoints,sqrt.(couplingRadius[:,k]),couplingRadius[:,k]), '\t')
        end
    end

end



function setup_system(
    ode, 
    system_name::AbstractString; 
    x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.0,0.0,0.0),
    attractor_fraction_coupling_list::Vector{Float64} = [0.01, 0.02, 0.05], 
    min_dist_tiling_list::Vector{Float64} = [0.25, 0.5], 
    attractor_fraction_tiling_list::Vector{Float64} = [0.001, 0.002]
    )
    
    setup_filename = string(system_name,".setup")

    #check if 
    if(isfile(DrWatson.datadir(system_name, setup_filename)))
        print("System already setup. Delete ",setup_filename," to repeat the setup.")
    else
        if(!isdir(DrWatson.datadir(system_name)))
            mkpath(DrWatson.datadir(system_name))
        end

        if(!isdir(DrWatson.plotsdir(system_name)))
            mkpath(DrWatson.plotsdir(system_name))
        end

        print("generating trajectories...\n")
        generate_attractor_plot(ode, system_name; x0)
        generate_attractor_sample(ode, system_name; x0)

        attractor_plot(system_name)

        print("generating attractor tilings...\n")
        generate_attractortiling_uniform(system_name, attractor_fraction_tiling_list = attractor_fraction_tiling_list)
        generate_attractortiling_equidistant(system_name, min_dist_tiling_list = min_dist_tiling_list)

        couplingregion_plot_all(system_name; attractor_fraction_tiling_list = attractor_fraction_tiling_list, min_dist_tiling_list = min_dist_tiling_list)

        print("generating coupling regions...\n")
        generate_couplingregions(system_name, attractor_fraction_coupling_list = attractor_fraction_coupling_list, min_dist_tiling_list = min_dist_tiling_list, attractor_fraction_tiling_list = attractor_fraction_tiling_list)

        touch(DrWatson.datadir(system_name,setup_filename))
    end

end


function simulate_nonlinear_trajectory(
    ode::Function, 
    system_name::AbstractString, 
    p::SystemParameters; 
    initial_state_index::Int64 = 1, 
    fixed_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.0,0.0,0.0), 
    delta_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.1,0.1,0.1), 
    simulation_time::Float64 = 2000.0, 
    transient_time::Float64 = 0.0, 
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4,
    #arguments for sampling and saving the trajectory
    sampling::Bool = false, 
    sampling_timestep::Float64 = 0.01, 
    output_filename::String = "nonlinear_trajectory.tsv"
    )   

    #set initial condition as first point from attractor sample
    x0_sample = (0.0,0.0,0.0)

    if initial_state_index > 0
        attractor_filename::String = string("attractor_sample.tsv")
        attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
        x0_sample = attractorpoints[initial_state_index,:]
    else
        x0_sample = fixed_x0
    end

    x0 = StaticArrays.SVector{6,Float64}(x0_sample[1], x0_sample[2], x0_sample[3], x0_sample[1]+delta_x0[1], x0_sample[2]+delta_x0[2], x0_sample[3]+delta_x0[3])

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            ode, 
            x0,
            p,
            diffeq = (
                save_everystep = false, 
                save_start = false,
                dtmax = max_integration_timestep,
                abstol = absolute_tolerance,
                reltol = relative_tolerance,
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)
    x0 = DynamicalSystems.current_state(ds)
    DynamicalSystems.reinit!(ds, x0, t0 = 0.0)

    if sampling
        u, t = DynamicalSystems.trajectory(ds, simulation_time; Δt = sampling_timestep)
        DelimitedFiles.writedlm(DrWatson.datadir(system_name, output_filename), hcat(t,u), '\t')
    else
        DynamicalSystems.step!(ds, simulation_time, true)
    end

    return (DynamicalSystems.current_time(ds), DynamicalSystems.current_state(ds))

end

function simulate_linearized_trajectory(
    ode::Function, 
    system_name::AbstractString, 
    p::SystemParameters; 
    initial_state_index::Int64 = 1, 
    fixed_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.0,0.0,0.0), 
    delta_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(1.0,1.0,1.0), 
    simulation_time::Float64 = 2000.0, 
    transient_time::Float64 = 600.0,  #2500.0
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4,
    #arguments for sampling and saving the trajectory
    sampling::Bool = false, 
    sampling_timestep::Float64 = 0.01, 
    output_filename::String = "linearized_system.tsv"
    )   

    #set initial condition as first point from attractor sample
    x0_sample = (0.0,0.0,0.0)

    if initial_state_index > 0
        attractor_filename::String = "attractor_sample.tsv"
        attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
        x0_sample = attractorpoints[initial_state_index,:]
    else
        x0_sample = fixed_x0
    end

    delta_x0_norm = sqrt( delta_x0[1]*delta_x0[1] + delta_x0[2]*delta_x0[2] + delta_x0[3]*delta_x0[3] )

    x0 = StaticArrays.SVector{7,Float64}(x0_sample[1], x0_sample[2], x0_sample[3], delta_x0[1]/delta_x0_norm, delta_x0[2]/delta_x0_norm, delta_x0[3]/delta_x0_norm, 0.0)   

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            ode, 
            x0,
            p,
            diffeq = (
                callback = linearizedODE_normalizecallback, 
                save_everystep = false, 
                save_start = false,
                dtmax = max_integration_timestep,
                abstol = absolute_tolerance,
                reltol = relative_tolerance
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)
    x_transient = DynamicalSystems.current_state(ds)
    x0_reinit = StaticArrays.setindex(x_transient, 0.0, 7) # Use setindex
    DynamicalSystems.reinit!(ds, x0_reinit, t0 = 0.0)

    if sampling
        u, t = DynamicalSystems.trajectory(ds, simulation_time; Δt = sampling_timestep)
        DelimitedFiles.writedlm(DrWatson.datadir(system_name, output_filename), hcat(t,u), '\t')
    else
        DynamicalSystems.step!(ds, simulation_time, true)
    end

    return (DynamicalSystems.current_time(ds), DynamicalSystems.current_state(ds))
    
end



function compute_maximum_lyapunov_exponent_linearized(
    ode::Function, 
    system_name::AbstractString, 
    p::SystemParameters; 
    x0::StaticArrays.SVector{7,Float64} = StaticArrays.SVector{7,Float64}(0.0,0.0,0.0,0.5,0.5,0.714,0.0), 
    simulation_time::Float64 = 100000.0, 
    transient_time::Float64 = 1000.0, 
    number_of_lyap_samples::Int64 = 10,
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4
    )   

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            ode, 
            x0,
            p,
            diffeq = (
                callback = linearizedODE_normalizecallback, 
                save_everystep = false, 
                save_start = false,
                dtmax = max_integration_timestep,
                abstol = absolute_tolerance,
                reltol = relative_tolerance
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)

    computed_lyap_exponents = StaticArrays.MVector{number_of_lyap_samples, Float64}(undef)

    for k in 1:number_of_lyap_samples 
        x_transient = DynamicalSystems.current_state(ds)
        
        x0_reinit = StaticArrays.setindex(x_transient, 0.0, 7) # Use setindex
        DynamicalSystems.reinit!(ds, x0_reinit, t0 = 0.0)

        DynamicalSystems.step!(ds, simulation_time)

        computed_lyap_exponents[k] = DynamicalSystems.current_state(ds)[7]/DynamicalSystems.current_time(ds)
    end

    #print("Estimated Lyapunov exponent: ", Statistics.mean(computed_lyap_exponents), " +/- ", Statistics.std(computed_lyap_exponents)/sqrt(number_of_lyap_samples))
    return (Statistics.mean(computed_lyap_exponents), Statistics.std(computed_lyap_exponents), Statistics.std(computed_lyap_exponents)/sqrt(number_of_lyap_samples))

end


function scan_maximum_lyapunov_exponent(
    ode::Function, 
    system_name::AbstractString, 
    p::SystemParameters; 
    output_filename_ovewrite::String = "", 
    coupling_strength_range::Vector{AbstractVector{Float64}} = [range(0.0, 5.0, step=0.01), zeros(501), zeros(501)], 
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

    lyap_scan_filename::String = string("lyap_scan.tsv")
    if(output_filename_ovewrite != "")
        lyap_scan_filename = output_filename_ovewrite
    end

    x0_sample = (0.0,0.0,0.0)

    if initial_state_index > 0
        attractor_filename::String = "attractor_sample.tsv"
        attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
        x0_sample = attractorpoints[initial_state_index,:]
    else
        x0_sample = fixed_x0
    end

    delta_x0_norm = sqrt( delta_x0[1]*delta_x0[1] + delta_x0[2]*delta_x0[2] + delta_x0[3]*delta_x0[3] )
    x0 = StaticArrays.SVector{7,Float64}(x0_sample[1], x0_sample[2], x0_sample[3], delta_x0[1]/delta_x0_norm, delta_x0[2]/delta_x0_norm, delta_x0[3]/delta_x0_norm, 0.0)   

    lyap_estimates::Matrix{Float64} = Matrix{Float64}(undef, size(coupling_strength_range[1],1), 6)

    if(show_progress)
        @time ProgressMeter.@showprogress Threads.@threads for k in 1:size(coupling_strength_range[1],1)
            alpha_x = coupling_strength_range[1][k]
            alpha_y = coupling_strength_range[2][k]
            alpha_z = coupling_strength_range[3][k]

            test_params = SystemParameters(
                StaticArrays.@SVector[alpha_x,alpha_y,alpha_z], 
                p.numberofcouplingregions, 
                p.couplingpoints, 
                p.couplingradii)

            mean, std, ste = compute_maximum_lyapunov_exponent_linearized(
                ode, system_name, 
                test_params; 
                x0 = x0, 
                simulation_time = simulation_time, 
                transient_time = transient_time, 
                number_of_lyap_samples = number_of_lyap_samples,
                max_integration_timestep = max_integration_timestep, 
                relative_tolerance = relative_tolerance,
                absolute_tolerance = absolute_tolerance
                )

            lyap_estimates[k,:] = vcat([alpha_x, alpha_y, alpha_z], [mean, std, ste])
        end
    else
        Threads.@threads for k in 1:size(coupling_strength_range[1],1)
            alpha_x = coupling_strength_range[1][k]
            alpha_y = coupling_strength_range[2][k]
            alpha_z = coupling_strength_range[3][k]
            
            test_params = SystemParameters(
                StaticArrays.@SVector[alpha_x,alpha_y,alpha_z], 
                p.numberofcouplingregions, 
                p.couplingpoints, 
                p.couplingradii)

            mean, std, ste = compute_maximum_lyapunov_exponent_linearized(
                ode, system_name, 
                test_params; 
                x0 = x0, 
                simulation_time = simulation_time, 
                transient_time = transient_time, 
                number_of_lyap_samples = number_of_lyap_samples,
                max_integration_timestep = max_integration_timestep, 
                relative_tolerance = relative_tolerance,
                absolute_tolerance = absolute_tolerance
                )

            lyap_estimates[k,:] = vcat([alpha], [mean, std, ste])
        end
    end

    DelimitedFiles.writedlm(DrWatson.datadir(system_name, lyap_scan_filename), lyap_estimates, '\t')

end

function compute_coupling_region_lyap_effect(
    ode::Function, 
    system_name::AbstractString, 
    coupling_region_filename::String, 
    coupling_strength::StaticArrays.SVector{3, Float64};
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

    couplingregion_data::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    couplingregion_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name, coupling_region_filename), '\t', Float64)

    number_of_coupling_regions::Int64 = number_of_coupling_region_samples
    if number_of_coupling_region_samples  == -1 || number_of_coupling_region_samples > size(couplingregion_data,1)
        number_of_coupling_regions = size(couplingregion_data,1)
    end

    x0_sample = (0.0,0.0,0.0)

    if initial_state_index > 0
        attractor_filename::String = "attractor_sample.tsv"
        attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
        x0_sample = attractorpoints[initial_state_index,:]
    else
        x0_sample = fixed_x0
    end

    delta_x0_norm = sqrt( delta_x0[1]*delta_x0[1] + delta_x0[2]*delta_x0[2] + delta_x0[3]*delta_x0[3] )
    x0 = StaticArrays.SVector{7,Float64}(x0_sample[1], x0_sample[2], x0_sample[3], delta_x0[1]/delta_x0_norm, delta_x0[2]/delta_x0_norm, delta_x0[3]/delta_x0_norm, 0.0)   

    uncoupled_params = SystemParameters(StaticArrays.@SVector[0.0, 0.0, 0.0], 0, StaticArrays.SVector{3, Float64}[], Float64[])
    base_mean, base_std, base_ste = compute_maximum_lyapunov_exponent_linearized(
        ode, 
        system_name, 
        uncoupled_params; 
        x0 = x0, 
        simulation_time = simulation_time, 
        transient_time = transient_time, 
        number_of_lyap_samples = number_of_lyap_samples,
        max_integration_timestep = max_integration_timestep, 
        relative_tolerance = relative_tolerance,
        absolute_tolerance = absolute_tolerance
        )

    couplingregion_lyap_estimates::Matrix{Float64} = Matrix{Float64}(undef, number_of_coupling_regions, 11)
    
    if(show_progress)
        @time ProgressMeter.@showprogress Threads.@threads for k in 1:number_of_coupling_regions
            test_region = StaticArrays.SVector{3, Float64}(couplingregion_data[k,1:3])
            test_params = SystemParameters(coupling_strength, 1, [test_region], [couplingregion_data[k,4]])
            mean, std, ste = compute_maximum_lyapunov_exponent_linearized(
                ode, 
                system_name, 
                test_params; 
                x0 = x0, 
                simulation_time = simulation_time, 
                transient_time = transient_time, 
                number_of_lyap_samples = number_of_lyap_samples,
                max_integration_timestep = max_integration_timestep, 
                relative_tolerance = relative_tolerance,
                absolute_tolerance = absolute_tolerance
                )

            couplingregion_lyap_estimates[k,:] = vcat(couplingregion_data[k,:], [mean, std, ste], [mean - base_mean, sqrt(std*std+base_std*base_std), ste+base_ste])
        end
    else
        Threads.@threads for k in 1:number_of_coupling_region_samples
            test_region = StaticArrays.SVector{3, Float64}(couplingregion_data[k,1:3])
            test_params = SystemParameters(coupling_strength, 1, [test_region], [couplingregion_data[k,4]])
            mean, std, ste = compute_maximum_lyapunov_exponent_linearized(
                ode, 
                system_name, 
                test_params; 
                x0 = x0, 
                simulation_time = simulation_time, 
                transient_time = transient_time, 
                number_of_lyap_samples = number_of_lyap_samples,
                max_integration_timestep = max_integration_timestep, 
                relative_tolerance = relative_tolerance,
                absolute_tolerance = absolute_tolerance
                )

            couplingregion_lyap_estimates[k,:] = vcat(couplingregion_data[k,:], [mean, std, ste], [mean - base_mean, sqrt(std*std+base_std*base_std), ste+base_ste])
        end
    end

    output_filename = string(coupling_region_filename[1:end-4], "__lyap_effect__alpha_",Printf.@sprintf("%.2f",coupling_strength[1]),"_",Printf.@sprintf("%.2f",coupling_strength[2]),"_",Printf.@sprintf("%.2f",coupling_strength[3]),".tsv")
    DelimitedFiles.writedlm(DrWatson.datadir(system_name, output_filename), couplingregion_lyap_estimates, '\t')
end

function find_next_optimal_coupling_region(
    ode::Function, 
    system_name::AbstractString, 
    p::SystemParameters, 
    coupling_region_filename::String; 
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

    couplingregion_data::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    couplingregion_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name, coupling_region_filename), '\t', Float64)

    number_of_coupling_regions = number_of_coupling_region_samples
    if number_of_coupling_region_samples  == -1 || number_of_coupling_region_samples > size(couplingregion_data,1)
        number_of_coupling_regions = size(couplingregion_data,1)
    end

    x0_sample = (0.0,0.0,0.0)

    if initial_state_index > 0
        attractor_filename::String = "attractor_sample.tsv"
        attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
        x0_sample = attractorpoints[initial_state_index,:]
    else
        x0_sample = fixed_x0
    end

    delta_x0_norm = sqrt( delta_x0[1]*delta_x0[1] + delta_x0[2]*delta_x0[2] + delta_x0[3]*delta_x0[3] )
    x0 = StaticArrays.SVector{7,Float64}(x0_sample[1], x0_sample[2], x0_sample[3], delta_x0[1]/delta_x0_norm, delta_x0[2]/delta_x0_norm, delta_x0[3]/delta_x0_norm, 0.0)   

    couplingregion_lyap_estimates::Matrix{Float64} = Matrix{Float64}(undef, number_of_coupling_regions, 8)
    
    if(show_progress)
        @time ProgressMeter.@showprogress Threads.@threads for k in 1:number_of_coupling_regions
            test_region = StaticArrays.SVector{3, Float64}(couplingregion_data[k,1:3])
            test_params = SystemParameters(p.couplingstrength, p.numberofcouplingregions + 1, vcat(p.couplingpoints,[test_region]), vcat(p.couplingradii,[couplingregion_data[k,4]]))
            mean, std, ste = compute_maximum_lyapunov_exponent_linearized(
                ode, system_name, 
                test_params; 
                x0 = x0, 
                simulation_time = simulation_time, 
                transient_time = transient_time, 
                number_of_lyap_samples = number_of_lyap_samples,
                max_integration_timestep = max_integration_timestep, 
                relative_tolerance = relative_tolerance,
                absolute_tolerance = absolute_tolerance
                )

            couplingregion_lyap_estimates[k,:] = vcat(couplingregion_data[k,:], [mean, std, ste])
            if couplingregion_data[k,1:3] in p.couplingpoints
                couplingregion_lyap_estimates[k,6] = Inf
            end
        end
    else
        Threads.@threads for k in 1:number_of_coupling_regions
            test_region = StaticArrays.SVector{3, Float64}(couplingregion_data[k,1:3])
            test_params = SystemParameters(p.couplingstrength, p.numberofcouplingregions + 1, vcat(p.couplingpoints,[test_region]), vcat(p.couplingradii,[couplingregion_data[k,4]]))
            mean, std, ste = compute_maximum_lyapunov_exponent_linearized(
                ode, system_name, 
                test_params; 
                x0 = x0, 
                simulation_time = simulation_time, 
                transient_time = transient_time, 
                number_of_lyap_samples = number_of_lyap_samples,
                max_integration_timestep = max_integration_timestep, 
                relative_tolerance = relative_tolerance,
                absolute_tolerance = absolute_tolerance
                )

            couplingregion_lyap_estimates[k,:] = vcat(couplingregion_data[k,:], [mean, std, ste])
            if couplingregion_data[k,1:3] in p.couplingpoints
                couplingregion_lyap_estimates[k,6] = Inf
            end
        end
    end

    best_coupling_region_index = argmin(couplingregion_lyap_estimates[:,6])
    best_coupling_region = couplingregion_lyap_estimates[best_coupling_region_index,:]

    return best_coupling_region
end


function simulate_nonlinear_optimal_coupling_trajectory(
    ode::Function, 
    system_name::AbstractString, 
    optimal_coupling_regions_filename::String, 
    coupling_strength::StaticArrays.SVector{3, Float64}; 
    number_of_coupling_regions::Int64 = 10, 
    initial_state_index::Int64 = 1, 
    fixed_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.0,0.0,0.0), 
    delta_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.1,0.1,0.1), 
    simulation_time::Float64 = 2000.0, 
    transient_time::Float64 = 0.0, 
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4,
    #arguments for sampling and saving the trajectory 
    sampling::Bool = false, 
    sampling_timestep::Float64 = 0.01, 
    output_suffix::String = ""
    )   
    
    optimal_coupling_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name, optimal_coupling_regions_filename), '\t', Float64)

    number_of_optimal_coupling_points = min(size(optimal_coupling_data,1), number_of_coupling_regions)
    optimal_coupling_points = [StaticArrays.SVector{3, Float64}(optimal_coupling_data[i, 1:3]) for i in 1:number_of_optimal_coupling_points]
    optimal_coupling_radii = optimal_coupling_data[:, 4]

    optimal_coupling_parameters = SystemParameters(coupling_strength, number_of_optimal_coupling_points, optimal_coupling_points, optimal_coupling_radii)

    output_filename = string("nonlinear__",optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_optimal_coupling_points,"__trajectory_",output_suffix,".tsv")

    simulate_nonlinear_trajectory(
        ode,
        system_name,
        optimal_coupling_parameters; 
        initial_state_index = initial_state_index,
        fixed_x0 = fixed_x0,
        delta_x0 = delta_x0,
        simulation_time = simulation_time,
        transient_time = transient_time,
        max_integration_timestep = max_integration_timestep,
        absolute_tolerance = absolute_tolerance,
        relative_tolerance = relative_tolerance,
        sampling = sampling, 
        sampling_timestep = sampling_timestep,
        output_filename = output_filename
        )
end

function simulate_linearized_optimal_coupling_trajectory(
    ode::Function, 
    system_name::AbstractString, 
    optimal_coupling_regions_filename::String, 
    coupling_strength::StaticArrays.SVector{3, Float64}; 
    number_of_coupling_regions::Int64 = 10, 
    initial_state_index::Int64 = 1, 
    fixed_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(0.0,0.0,0.0), 
    delta_x0::StaticArrays.SVector{3,Float64} = StaticArrays.SVector{3,Float64}(1.0, 1.0, 1.0), 
    simulation_time::Float64 = 2000.0, 
    transient_time::Float64 = 2500.0, 
    max_integration_timestep::Float64 = 0.01, 
    absolute_tolerance::Float64 = 1e-6, 
    relative_tolerance::Float64 = 1e-4,
    #arguments for sampling and saving the trajectory 
    sampling::Bool = false, 
    sampling_timestep::Float64 = 0.01, 
    output_suffix::String = ""
    )   
    
    optimal_coupling_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name, optimal_coupling_regions_filename), '\t', Float64)

    number_of_optimal_coupling_points = min(size(optimal_coupling_data,1), number_of_coupling_regions)
    optimal_coupling_points = [StaticArrays.SVector{3, Float64}(optimal_coupling_data[i, 1:3]) for i in 1:number_of_optimal_coupling_points]
    optimal_coupling_radii = optimal_coupling_data[:, 4]

    optimal_coupling_parameters = SystemParameters(coupling_strength, number_of_optimal_coupling_points, optimal_coupling_points, optimal_coupling_radii)

    output_filename = string("linearized__",optimal_coupling_regions_filename[1:end-4],"_ncr_",number_of_optimal_coupling_points,"__trajectory_",output_suffix,".tsv")

    simulate_linearized_trajectory(
        ode,
        system_name,
        optimal_coupling_parameters; 
        initial_state_index = initial_state_index,
        fixed_x0 = fixed_x0,
        delta_x0 = delta_x0,
        simulation_time = simulation_time,
        transient_time = transient_time,
        max_integration_timestep = max_integration_timestep,
        absolute_tolerance = absolute_tolerance,
        relative_tolerance = relative_tolerance,
        sampling = sampling, 
        sampling_timestep = sampling_timestep,
        output_filename = output_filename
    )
end


function compute_coupling_fraction( 
    system_name::AbstractString, 
    optimal_coupling_regions_filename::String
    )   

    couplingregion_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name, optimal_coupling_regions_filename), '\t', Float64)
    total_coupling_regions = size(couplingregion_data,1)

    attractor_filename::String = "attractor_sample.tsv"
    attractorpoints = DelimitedFiles.readdlm(DrWatson.datadir(system_name,attractor_filename), '\t', Float64)
    total_attractor_points = size(attractorpoints,1)

    coupling_fractions::Matrix{Float64} = hcat(couplingregion_data,zeros(total_coupling_regions))

    for number_of_coupling_regions::Int64 in 1:total_coupling_regions
        attractor_points_inside_region = 0
        for k::Int64 in 1:total_attractor_points
            for j::Int64 in 1:number_of_coupling_regions
                if( (attractorpoints[k,1]-couplingregion_data[j,1])^2 + (attractorpoints[k,2]-couplingregion_data[j,2])^2 + (attractorpoints[k,3]-couplingregion_data[j,3])^2 < couplingregion_data[j,4]^2 )
                    attractor_points_inside_region += 1
                    break
                end
            end
        end
        
        coupling_fractions[number_of_coupling_regions,9] = attractor_points_inside_region / total_attractor_points
    end

    DelimitedFiles.writedlm(DrWatson.datadir(system_name, optimal_coupling_regions_filename), coupling_fractions, '\t')

end


function is_in_coupling_region(x::Vector{Float64}, coupling_region_filename::String, numberofcouplingregions::Int64, )

    coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler",coupling_region_filename), '\t', Float64)

    for k::Int in 1:numberofcouplingregions
        if( (x[1]-coupling_regions_data[k,1])^2 + (x[2]-coupling_regions_data[k,2])^2 + (x[3]-coupling_regions_data[k,3])^2 < coupling_regions_data[k,4]^2 )
            return true
        end
    end

    return false
end