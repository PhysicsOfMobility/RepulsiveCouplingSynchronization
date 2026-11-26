using DrWatson
@quickactivate :ChaosSynchronization

import DelimitedFiles
import StaticArrays
import LinearAlgebra

using CairoMakie


#define TU Dresden Colors
function RGBA(r::Int64,g::Int64,b::Int64,a::Int64)
    RGBAf(r/255.0, g/255.0, b/255.0, a/255.0)
end
function RGB(r::Int64,g::Int64,b::Int64)
    RGBA(r, g, b, 255)
end

TUD_darkblue = RGB(0,20,80)
TUD_deepblue = RGB(0,0,140)
TUD_darkgreen = RGB(0,60,15)
TUD_deepgreen = RGB(0,110,0)
TUD_darkred = RGB(80,0,20)
TUD_deepred = RGB(140,0,0)


TUD_lightblue1 = RGB(47,87,178)
TUD_lightblue2 = RGB(170,190,220)
TUD_purple1 = RGB(115,105,190)
TUD_purple2 = RGB(200,200,255)
TUD_magenta1 = RGB(188,21,137)
TUD_magenta2 = RGB(255,185,255)
TUD_red1 = RGB(210,12,65)
TUD_red2 = RGB(255,170,165)
TUD_orange1 = RGB(200,80,0)
TUD_orange2 = RGB(255,190,120)
TUD_yellow1 = RGB(255,199,0)
TUD_yellow2 = RGB(255,228,131)
TUD_oliv1 = RGB(118,122,35)
TUD_oliv2 = RGB(210,220,70)
TUD_green1 = RGB(0,125,75)
TUD_green2 = RGB(140,230,170)
TUD_turquoise1 = RGB(10,119,127)
TUD_turquoise2 = RGB(140,230,215)

TUD_black = RGB(0,0,0)
TUD_white = RGB(255,255,255)
TUD_grey100 = RGB(50,63,75)
TUD_grey80 = RGB(86,99,113)
TUD_grey60 = RGB(125,136,148)
TUD_grey40 = RGB(165,174,184)
TUD_grey20 = RGB(208,213,220)
TUD_grey10 = RGB(231,233,237)


unit_per_mm = 72.0/25.4 #72 pt per inch / 25.4 mm per inch
#unit_per_mm = 90.0/25.4 
mm_per_unit = 1.0/unit_per_mm

#sizes for full column figures (183mm from Nature)
TUD_theme = Theme(
    pt_per_unit = 1,
    px_per_unit = 2.666,
    #size=(530,258),    #187mm x 91mm in pt mit 72 pt per inch--> 187/25.4*72,91/25.4*72 for closest match
    #size=(530,300),    #higher in y-dimension to better fit picture aspect ratios
    #size=(663,375),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
    size=(660,660),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
    figure_padding=1*unit_per_mm,  #1mm padding
    fontsize = 10,  #in pt
    fonts = (;
        regular = "Noto Sans", 
        bold = "Noto Sans Bold", 
        italic = "Noto Sans Italic", 
        bold_italic = "Noto Sans Bold Italic", 
        tex_text_regular = "Computer Modern",
        tex_text_bold = "Computer Modern Bold",
        tex_text_italic = "Computer Modern Italic",
        tex_text_bold_italic = "Computer Modern Bold Italic",
        tex_math = "Computer Modern Italic",
        tex_math_bold = "Computer Modern Bold Italic"
    ),  
    backgroundcolor=TUD_white,
    xgridvisible = false,
    ygridvisible = false
#    palette = (
#        color = (   
#            TUD_darkblue, TUD_darkgreen, TUD_darkred, 
#            TUD_lightblue1, TUD_green1, TUD_red1, 
#            TUD_purple1, TUD_turquoise1, TUD_orange1,
#            TUD_magenta1, TUD_oliv1, TUD_yellow1,
#            TUD_deepblue, TUD_deepgreen, TUD_deepred, 
#            TUD_lightblue2, TUD_green2, TUD_red2, 
#            TUD_purple2, TUD_turquoise2, TUD_orange2,
#            TUD_magenta2, TUD_oliv2, TUD_yellow2,
#        ),
#    )
)

set_theme!(TUD_theme)




#set up all data required for later figures


coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","couplingregions_equidist__dist_0.50__cfrac_0.01.tsv"), '\t', Float64)
optimal_coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00.tsv"), '\t', Float64)


coupling_region_example_1a_index = 1243 
coupling_region_example_2a_index = 556

coupling_region_example_1b_index = 1243 
coupling_region_example_2b_index = 564

default_color = TUD_darkblue
uncoupled_color = TUD_grey60
optimal_coupling_region_color = TUD_green1
optimal_coupling_region_timeseries_color = RGB(70,178,125)
coupling_region_example_1_color = TUD_orange1
coupling_region_example_2_color = TUD_orange2
coupling_region_example_3_color = TUD_red1

attractor_plot_range = 1:50000
number_of_optimal_coupling_regions = 15


synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
init_state = synch_error_data_opt_full[1,2:8]

systemparameters_uncoupled = SystemParameters(
    StaticArrays.@SVector[0.0, 0.0, 0.0], 
    0, 
    StaticArrays.SVector{3, Float64}[], 
    Float64[]
)

simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4], init_state[5], init_state[6]], transient_time = 0.0, sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 0.01)


#
#synch_error_data_opt_full_nonlinear = DelimitedFiles.readdlm(DrWatson.datadir("roessler","nonlinear__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
#init_state_nonlinear  = synch_error_data_opt_full[1,2:7]
#
#simulate_nonlinear_trajectory(roesslerODE_coupled,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4]-init_state[1], init_state[5]-init_state[2], init_state[6]-init_state[3]], simulation_time = 1.0e6, transient_time = 0.0, sampling = true, output_filename = "nonlinear__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 1.0)
#


#structure of figure

fig = Figure(
            #image size in Makie units converted with 300dpi from 183mm x 90mm (using Makies 0.75 pt_per_unit)
            #size = (ceil(183*25.4*300/0.75), ceil(90*25.4*300/0.75))
            figure_padding = (-150, -210, -115, -180)
            #backgroundcolor = :transparent
        )

#overall figure layout
group_fig = fig[1,1] = GridLayout()

panel_a = group_fig[1,1] = GridLayout()

#panel a: attractor plot optimal coupling regions

axis_panel_a = 
Axis3(
    panel_a[1, 1],
    aspect = :equal,
    azimuth = 1.35 * pi,
    elevation = 0.15 * pi,
    limits = (-10,14.9,-12,8,-2,25)
)

hidedecorations!(axis_panel_a)
hidespines!(axis_panel_a)

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_a, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
arrow_plot_range = 1:150:200000
arrows3d!(
    axis_panel_a,
    synch_error_data_opt_full[arrow_plot_range,2], 
    synch_error_data_opt_full[arrow_plot_range,3], 
    synch_error_data_opt_full[arrow_plot_range,4], 
    synch_error_data_opt_full[arrow_plot_range,5],
    synch_error_data_opt_full[arrow_plot_range,6],
    synch_error_data_opt_full[arrow_plot_range,7],
    shaftcolor = default_color, 
    tipcolor = default_color,
    alpha = 0.5,
    tiplength = 0.0,
    shaftradius = 0.05,
    normalize = true,
    lengthscale = 0.75,
    align = :center,
    transparency = true,
    rasterize = 10
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_a, 
    optimal_coupling_regions_data[1:number_of_optimal_coupling_regions,1],
    optimal_coupling_regions_data[1:number_of_optimal_coupling_regions,2],
    optimal_coupling_regions_data[1:number_of_optimal_coupling_regions,3],
    color = optimal_coupling_region_timeseries_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data[k,1], optimal_coupling_regions_data[k,2], optimal_coupling_regions_data[k,3]) for k in 1:number_of_optimal_coupling_regions]
meshscatter!(axis_panel_a,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data[1:number_of_optimal_coupling_regions,4], 
    color = optimal_coupling_region_timeseries_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)



save(plotsdir("FIG_title.pdf"), fig)
save(plotsdir("FIG_title.png"), fig)


fig
