import DrWatson
import DelimitedFiles 
import Printf
import LinearAlgebra

using CairoMakie
CairoMakie.activate!()

function attractor_plot(system_name::AbstractString; kwargs...)

    attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,"attractor_plot.tsv"), '\t', Float64)
    attractor_sample_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,"attractor_sample.tsv"), '\t', Float64)

    # --- Plotting in 3D ---
    attractor_plot_fig = Figure(size = (1600, 1200))

    # Use Axis3 for a 3D plot
    ax = Axis3(attractor_plot_fig[1, 1]; 
        xlabel = "x", 
        ylabel = "y", 
        zlabel = "z",
        aspect = :equal,
        azimuth = 1.4 * pi,
        elevation = 0.15 * pi
    )

    lines!(ax, attractor_plot_data[:,2:4],
        linewidth = 0.25, 
        #color = RGBf(0.0/255, 0.0/255, 178.0/255),
        color = RGBf(0.0/255, 20.0/255, 140.0/255),
        #color = RGBf(47.0/255, 87.0/255, 178.0/255),
        alpha = 0.5
    )

    scatter!(
        ax, 
        attractor_sample_data,
        color = RGBf(47.0/255, 178.0/255, 87.0/255), 
        markersize = 2,      # Set a clear marker size
        strokewidth = 0.0,     # Remove stroke for cleaner appearance
        alpha = 0.1,
        transparency=true
    )

    #Save the figure
    save(DrWatson.plotsdir(system_name,"attractor.png"), attractor_plot_fig)
end



function couplingregion_plot(system_name::AbstractString; tiling_scheme::AbstractString = "equidistant", tiling_param::Float64 = 0.50)

    couplingregion_data::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    output_filename = ""

    attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,"attractor_plot.tsv"), '\t', Float64)

    if tiling_scheme == "equidist"
        tiling_dist = tiling_param
        tiling_filename = string("attractortiling_",tiling_scheme,"__dist_",Printf.@sprintf("%.2f",tiling_dist),".tsv")
        output_filename = string("attractortiling_",tiling_scheme,"__dist_",Printf.@sprintf("%.2f",tiling_dist),".png")

        couplingregion_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,tiling_filename), '\t', Float64)
    elseif tiling_scheme == "uniform"
        tiling_frac = tiling_param
        tiling_filename = string("attractortiling_",tiling_scheme,"__frac_",Printf.@sprintf("%.3f",tiling_frac),".tsv")
        output_filename = string("attractortiling_",tiling_scheme,"__frac_",Printf.@sprintf("%.3f",tiling_frac),".png")

        couplingregion_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,tiling_filename), '\t', Float64)
    else
        return
    end

    # --- Plotting in 3D ---
    couplingregion_plot_fig = Figure(size = (1600, 1200))

    ax = Axis3(couplingregion_plot_fig[1, 1]; 
        xlabel = "x", 
        ylabel = "y", 
        zlabel = "z",
        aspect = :equal,
        azimuth = 1.4 * pi,
        elevation = 0.15 * pi
    )

    lines!(ax, attractor_plot_data[:,2:4],
        linewidth = 0.25, 
        color = RGBf(47.0/255, 87.0/255, 178.0/255),
        alpha = 0.5,
        transparency = true
    )

    scatter!(
        ax, 
        couplingregion_data,
        color = RGBf(47.0/255, 178.0/255, 87.0/255), 
        markersize = 10,      # Set a clear marker size
        strokewidth = 0.0,     # Remove stroke for cleaner appearance
        alpha = 1.0,
        transparency=true
    )

    #Save the figure
    save(DrWatson.plotsdir(system_name, output_filename), couplingregion_plot_fig)
end


function couplingregion_plot_all(system_name::AbstractString; min_dist_tiling_list::Vector{Float64} = [0.5], attractor_fraction_tiling_list::Vector{Float64} = [0.002])
    for tiling_dist in min_dist_tiling_list
        couplingregion_plot(system_name; tiling_scheme = "equidist", tiling_param = tiling_dist)
    end

    for tiling_fraction in attractor_fraction_tiling_list
        couplingregion_plot(system_name; tiling_scheme = "uniform", tiling_param = tiling_fraction) 
    end
end


function optimal_coupling_regions_plot(system_name::AbstractString, optimal_coupling_regions_filename::String; number_of_coupling_regions::Int64 = 10)

    attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,"attractor_plot.tsv"), '\t', Float64)

    optimal_coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,optimal_coupling_regions_filename), '\t', Float64)
    number_of_coupling_regions_to_plot = min(size(optimal_coupling_regions_data,1),number_of_coupling_regions)

    # --- Plotting in 3D ---
    optimal_coupling_regions_plot_fig = Figure(size = (1600, 1200))

    ax = Axis3(optimal_coupling_regions_plot_fig[1, 1]; 
        xlabel = "x", 
        ylabel = "y", 
        zlabel = "z",
        aspect = :equal,
        azimuth = 1.4 * pi,
        elevation = 0.15 * pi
    )

    lines!(ax, attractor_plot_data[:,2:4],
        linewidth = 0.25, 
        color = RGBf(47.0/255, 87.0/255, 178.0/255),
        alpha = 0.5,
        transparency = true
    )

    scatter!(
        ax, 
        optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,1],
        optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,2],
        optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,3],
        color = RGBf(47.0/255, 178.0/255, 87.0/255), 
        markersize = 10,      # Set a clear marker size
        strokewidth = 0.0,     # Remove stroke for cleaner appearance
        alpha = 1.0,
        transparency=true
    )

    spherecenters = [Point(optimal_coupling_regions_data[k,1], optimal_coupling_regions_data[k,2], optimal_coupling_regions_data[k,3]) for k in 1:(number_of_coupling_regions_to_plot)]

    meshscatter!(ax,
        spherecenters,
        marker = Sphere(Point(0.0,0.0,0.0),1),
        markersize = optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,4], 
        color = RGBf(47.0/255, 178.0/255, 87.0/255),
        alpha = 0.5,
        transparency=true
    )

    #Save the figure
    output_filename = string(optimal_coupling_regions_filename[1:end-4], "_ncr_",number_of_coupling_regions,".png")
    save(DrWatson.plotsdir(system_name, output_filename), optimal_coupling_regions_plot_fig)
end



function attractor_difference_direction_plot(system_name::AbstractString; difference_filename_ovewrite::String = "", output_filename_overwrite::String = "")

    attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,"attractor_plot.tsv"), '\t', Float64)
    difference_filename::String = "linearized__uncoupled_trajectory.tsv"
    
    if difference_filename_ovewrite != ""
        difference_filename = difference_filename_ovewrite
    end

    attractor_difference_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,difference_filename), '\t', Float64)

    # --- Plotting in 3D ---
    attractor_plot_fig = Figure(size = (1600, 1200))

    # Use Axis3 for a 3D plot
    ax = Axis3(attractor_plot_fig[1, 1]; 
        xlabel = "x", 
        ylabel = "y", 
        zlabel = "z",
        aspect = :equal,
        azimuth = 1.4 * pi,
        elevation = 0.15 * pi
    )

    lines!(ax, attractor_plot_data[:,2:4],
        linewidth = 0.25, 
        #color = RGBf(0.0/255, 0.0/255, 178.0/255),
        color = RGBf(0.0/255, 20.0/255, 140.0/255),
        #color = RGBf(47.0/255, 87.0/255, 178.0/255),
        alpha = 0.5,
        transparency = true
    )

    arrows3d!(
        ax,
        attractor_difference_data[:,2], 
        attractor_difference_data[:,3], 
        attractor_difference_data[:,4], 
        attractor_difference_data[:,5],
        attractor_difference_data[:,6],
        attractor_difference_data[:,7],
        shaftcolor = RGBf(47.0/255, 178.0/255, 87.0/255), 
        tipcolor = RGBf(47.0/255, 178.0/255, 87.0/255),
        tiplength = 0.0,
        shaftradius = 0.015,
        normalize = true,
        lengthscale = 0.25,
        align = :center,
        transparency = true
    )

    #Save the figure
    output_filename = ""
    if output_filename_overwrite != ""
        output_filename = output_filename_overwrite
    else
        output_filename = string(difference_filename[1:end-4], "__difference_vectors.png")
    end
    save(DrWatson.plotsdir(system_name,output_filename), attractor_plot_fig)
end


function attractor_difference_direction_with_optimal_coupling_regions_plot(system_name::AbstractString, optimal_coupling_regions_filename::String, difference_trajectory_filename::String = ""; number_of_coupling_regions::Int64 = 10, )

    attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,"attractor_plot.tsv"), '\t', Float64)

    attractor_difference_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,difference_trajectory_filename), '\t', Float64)

    # --- Plotting in 3D ---
    attractor_plot_fig = Figure(size = (1600, 1200))

    # Use Axis3 for a 3D plot
    ax = Axis3(attractor_plot_fig[1, 1]; 
        xlabel = "x", 
        ylabel = "y", 
        zlabel = "z",
        aspect = :equal,
        azimuth = 1.4 * pi,
        elevation = 0.15 * pi
    )

    lines!(ax, attractor_plot_data[:,2:4],
        linewidth = 0.25, 
        #color = RGBf(0.0/255, 0.0/255, 178.0/255),
        color = RGBf(0.0/255, 20.0/255, 140.0/255),
        #color = RGBf(47.0/255, 87.0/255, 178.0/255),
        alpha = 0.5,
        transparency = true
    )

    arrows3d!(
        ax,
        attractor_difference_data[:,2], 
        attractor_difference_data[:,3], 
        attractor_difference_data[:,4], 
        attractor_difference_data[:,5],
        attractor_difference_data[:,6],
        attractor_difference_data[:,7],
        shaftcolor = RGBf(47.0/255, 178.0/255, 87.0/255), 
        tipcolor = RGBf(47.0/255, 178.0/255, 87.0/255),
        tiplength = 0.0,
        shaftradius = 0.015,
        normalize = true,
        lengthscale = 0.25,
        align = :center,
        transparency = true
    )


    optimal_coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,optimal_coupling_regions_filename), '\t', Float64)
    number_of_coupling_regions_to_plot = min(size(optimal_coupling_regions_data,1),number_of_coupling_regions)

    scatter!(
        ax, 
        optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,1],
        optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,2],
        optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,3],
        color = RGBf(47.0/255, 178.0/255, 87.0/255), 
        markersize = 10,      # Set a clear marker size
        strokewidth = 0.0,     # Remove stroke for cleaner appearance
        alpha = 1.0,
        transparency=true
    )

    spherecenters = [Point(optimal_coupling_regions_data[k,1], optimal_coupling_regions_data[k,2], optimal_coupling_regions_data[k,3]) for k in 1:number_of_coupling_regions_to_plot]

    meshscatter!(ax,
        spherecenters,
        marker = Sphere(Point(0.0,0.0,0.0),1),
        markersize = optimal_coupling_regions_data[1:number_of_coupling_regions_to_plot,4], 
        #color = RGBf(47.0/255, 178.0/255, 87.0/255),
        color = RGBf(47.0/255, 178.0/255, 255.0/255),
        alpha = 0.5,
        transparency=true
    )

    #Save the figure
    output_filename = string(difference_trajectory_filename[1:end-4], "__difference_vectors_opt.png")
    save(DrWatson.plotsdir(system_name,output_filename), attractor_plot_fig)
end


function couplingregion_lyap_effect_plot(system_name::AbstractString, coupling_region_lyap_effect_filename::String)

    output_filename = ""

    attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,"attractor_plot.tsv"), '\t', Float64)

    couplingregion_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,coupling_region_lyap_effect_filename), '\t', Float64)

    # --- Plotting in 3D ---
    couplingregion_lyap_effect_plot_fig = Figure(size = (1600, 1200))

    ax = Axis3(couplingregion_lyap_effect_plot_fig[1, 1]; 
        xlabel = "x", 
        ylabel = "y", 
        zlabel = "z",
        aspect = :equal,
        azimuth = 1.4 * pi,
        elevation = 0.15 * pi
    )

    lines!(ax, attractor_plot_data[:,2:4],
        linewidth = 0.25, 
        color = RGBf(47.0/255, 87.0/255, 178.0/255),
        alpha = 0.5,
        transparency = true
    )

    scatter!(
        ax, 
        couplingregion_data[:,1],
        couplingregion_data[:,2],
        couplingregion_data[:,3],
        color = couplingregion_data[:,9], 
        colormap = :berlin,
        colorscale = Makie.Symlog10(0.001),
        colorrange = (-0.25, 0.25),
        markersize = 30 * 2.5e-4 ./ couplingregion_data[:,11],      # Set a clear marker size
        strokewidth = 0.0,     
        alpha = 1.0,
        transparency=true
    )

    #Save the figure
    output_filename = string(coupling_region_lyap_effect_filename[1:end-4], ".png")
    save(DrWatson.plotsdir(system_name, output_filename), couplingregion_lyap_effect_plot_fig)
end




function synch_error_plot_nonlinear(system_name::AbstractString, trajectory_filename::String; kwargs...)

    trajectory_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,trajectory_filename), '\t', Float64)

    # --- Plotting in 3D ---
    synch_error_plot_fig = Figure(size = (1600, 1200))

    ax = Axis(synch_error_plot_fig[1, 1]; 
        xlabel = "time", 
        ylabel = "synchronization error", 
        aspect = 2,
        yscale = log10,
        xlabelsize = 28,
        ylabelsize = 28,
        xticklabelsize = 20,
        yticklabelsize = 20,
        limits = (0.0, 100.0, 1e-6, 1e3) # 'nothing' lets y-axis auto-scale
    )

    lines!(ax, 
        trajectory_data[:,1],
        LinearAlgebra.norm.(eachrow(trajectory_data[:,5:7] - trajectory_data[:,2:4])),
        linewidth = 2, 
        color = RGBf(0.0/255, 20.0/255, 140.0/255)
    )

    #Save the figure
    output_filename = string(trajectory_filename[1:end-4], "__syncherror.png")
    save(DrWatson.plotsdir(system_name, output_filename), synch_error_plot_fig)
end

function synch_error_plot_linearized(system_name::AbstractString, trajectory_filename::String; kwargs...)

    trajectory_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,trajectory_filename), '\t', Float64)

    # --- Plotting in 3D ---
    synch_error_plot_fig = Figure(size = (1600, 1200))

    ax = Axis(synch_error_plot_fig[1, 1]; 
        xlabel = "time", 
        ylabel = "synchronization error", 
        aspect = 2,
        yscale = log10,
        xlabelsize = 28,
        ylabelsize = 28,
        xticklabelsize = 20,
        yticklabelsize = 20,
        limits = (0.0, 100.0, 1e-6, 1e3) # 'nothing' lets y-axis auto-scale
    )

    lines!(ax, 
        trajectory_data[:,1],
        exp.(trajectory_data[:,8]),
        linewidth = 2, 
        color = RGBf(0.0/255, 20.0/255, 140.0/255)
    )

    #Save the figure
    output_filename = string(trajectory_filename[1:end-4], "__syncherror.png")
    save(DrWatson.plotsdir(system_name, output_filename), synch_error_plot_fig)
end


function lyap_scan_plot(system_name::AbstractString, lyap_scan_filename::String; scan_dimension::Int64 = 1, kwargs...)

    lyap_scan_data = DelimitedFiles.readdlm(DrWatson.datadir(system_name,lyap_scan_filename), '\t', Float64)

    # --- Plotting in 3D ---
    lyap_scan_plot_fig = Figure(size = (1600, 1200))

    ax = Axis(lyap_scan_plot_fig[1, 1]; 
        xlabel = "coupling strength", 
        ylabel = "lyapunov exponent", 
        aspect = 2,
        xlabelsize = 28,
        ylabelsize = 28,
        xticklabelsize = 20,
        yticklabelsize = 20,
        limits = (0.0, 5.0, nothing, nothing) # 'nothing' lets y-axis auto-scale
    )

    lines!(ax, 
        lyap_scan_data[:,scan_dimension],
        lyap_scan_data[:,4],
        linewidth = 2, 
        color = RGBf(0.0/255, 20.0/255, 140.0/255)
    )

    #Save the figure
    output_filename = string(lyap_scan_filename[1:end-4], ".png")
    save(DrWatson.plotsdir(system_name, output_filename), lyap_scan_plot_fig)
end