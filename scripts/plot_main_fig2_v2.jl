using DrWatson
@quickactivate :ChaosSynchronization

import DelimitedFiles
import StaticArrays

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
    size=(663,300),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
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

#simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4], init_state[5], init_state[6]], transient_time = 0.0, sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 0.01)
simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4], init_state[5], init_state[6]], transient_time = 0.0, sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 0.05)


#
#synch_error_data_opt_full_nonlinear = DelimitedFiles.readdlm(DrWatson.datadir("roessler","nonlinear__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_10__trajectory_plot.tsv"), '\t', Float64)
#init_state_nonlinear  = synch_error_data_opt_full[1,2:7]
#
#simulate_nonlinear_trajectory(roesslerODE_coupled,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4]-init_state[1], init_state[5]-init_state[2], init_state[6]-init_state[3]], simulation_time = 1.0e6, transient_time = 0.0, sampling = true, output_filename = "nonlinear__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 1.0)
#


#structure of figure

fig = Figure(
            #image size in Makie units converted with 300dpi from 183mm x 90mm (using Makies 0.75 pt_per_unit)
            #size = (ceil(183*25.4*300/0.75), ceil(90*25.4*300/0.75))
        )

#overall figure layout
group_fig = fig[1,1] = GridLayout()

panel_a = group_fig[1,1] = GridLayout()
group_bcdefg = group_fig[1,2] = GridLayout()
panel_b = group_bcdefg[1,1] = GridLayout()
panel_c = group_bcdefg[2,1] = GridLayout()
panel_d = group_bcdefg[3,1] = GridLayout()
panel_e = group_bcdefg[1,2] = GridLayout()
panel_f = group_bcdefg[2,2] = GridLayout()
panel_g = group_bcdefg[3,2] = GridLayout()



colsize!(group_fig,1,Auto(0.425))
colsize!(group_fig,2,Auto(0.575))

colsize!(group_bcdefg,1,Auto(0.5))
colsize!(group_bcdefg,2,Auto(0.5))

rowsize!(group_bcdefg,1,Auto(0.33))
rowsize!(group_bcdefg,2,Auto(0.33))
rowsize!(group_bcdefg,3,Auto(0.33))

colgap!(group_fig, 4*unit_per_mm)
colgap!(group_bcdefg, 3*unit_per_mm)
rowgap!(group_bcdefg, 3*unit_per_mm)



#individual panels



#panel a: attractor plot optimal coupling regions

axis_panel_a = 
Axis3(
    panel_a[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.35 * pi,
    elevation = 0.15 * pi,
    limits = (-10,14.9,-12,8,-2,25),
    xticks = [-10,-5,0,5,10],
    xtickformat = values -> [L"-10",L"-5",L"0",L"5",L"10"],
    yticks = [-10,-5,0,5],
    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
    zticks = [0,5,10,15,20,25],
    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20",L"25"],
    protrusions = (25,0,10,0),
    xticklabelpad = 1,
    yticklabelpad = 1,
    zticklabelpad = 2,
    xlabeloffset = 15,
    ylabeloffset = 20,
    zlabeloffset = 20
    #ylabelpadding = -0.5,
    #xticklabelpad = -0.6,
    #yticklabelpad = 1,
    #yticks = [-0.5,-0.25,0],
    #ytickformat = values -> ["-0.50","-0.25","0"],
    #xtickalign = 1,
    #ytickalign = 1,    
    #ztickalign = 1
    #xgridvisible = false,
    #ygridvisible = false
)

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_a, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
arrow_plot_range = 1:15:20000
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
    color = optimal_coupling_region_color, 
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
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)


text!(axis_panel_a, 
    (0.085, 0.98), 
    text = L"\alpha = -1.0", 
    space = :relative, 
    align = (:left, :top), 
    color = :black, 
    fontsize = 14
    )    

text!(axis_panel_a, 
    (0.975, 0.98), 
    text = L"m = 15", 
    space = :relative, 
    align = (:right, :top), 
    color = :black,
    fontsize = 14
    )  



text!(axis_panel_a, 
    (0.95, 0.55), 
    text = L"\textbf{A}", 
    space = :relative,
    align = (:right, :top),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_a,
    (0.9, 0.5),
    (-0.1, -0.04),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    space = :relative,
)

text!(axis_panel_a, 
    (0.3, 0.56), 
    text = L"\textbf{B}", 
    space = :relative,
    align = (:right, :top),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_a,
    (0.28, 0.50),
    (0.00, -0.1),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    space = :relative,
)


text!(axis_panel_a, 
    (0.4, 0.14), 
    text = L"\textbf{C}", 
    space = :relative,
    align = (:right, :top),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_a,
    (0.425, 0.135),
    (0.1, 0.065),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    space = :relative,
)

Label(
    panel_a[1, 1, TopLeft()], 
    "a",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)



#panel b: optimal synch error

axis_panel_b = 
Axis(
    panel_b[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,5e3 ,1e-50, 1e10),
    #xticks = [0,200,400,600,800],
    #xtickformat = values -> [L"0","",L"400","",L"800",""],
    xticks = [0,1000,2000,3000,4000,5000],
    xtickformat = values -> [L"0",L"1000",L"2000",L"3000",L"4000",L"5000"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    #yticks = [1e-18,1e-12,1e-6,1e0,1e6,1e12,1e18,1e24],
    #ytickformat = values -> ["",L"10^{-12}","",L"10^0","",L"10^{12}","",L"10^{24}"],
    yticks = [1e-40,1e-20,1e0],
    ytickformat = values -> [L"10^{-40}",L"10^{-20}",L"10^0"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__uncoupled_trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_b,
    synch_error_data_uncoupled[:,1],exp.(synch_error_data_uncoupled[:,8]),
    color = uncoupled_color,
    linewidth = 1.5
)

synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_b,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_opt_full[:,1],exp.(synch_error_data_opt_full[:,8]),
    color = TUD_darkblue,
    linewidth = 1.5
)

lines!(
    axis_panel_b,
    [500,4500],[1e-13,1.3525e-34],
    color = TUD_darkblue,
    linewidth = 1
)

text!(
    axis_panel_b,
    (1500,1e-30),
    text = L"\lambda \approx -0.012",  #-0.012013082158072463
    rotation = -0.065*pi,
    color = TUD_darkblue
)

lines!(
    axis_panel_b,
    [0,1e4],[1e0,1e0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

Label(
    panel_b[1, 1, TopLeft()], 
    "b",
    fontsize = 12,
    font = :bold,
    padding = (-8*unit_per_mm, 10*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)





synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
synch_error_data_opt_full_active_coupling = [ 
    is_in_coupling_region(synch_error_data_opt_full[k,2:4], "couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00.tsv", number_of_optimal_coupling_regions) 
    for k in 1:10000
]


t_full = synch_error_data_opt_full[1:10000,1]
deltax_full = exp.(synch_error_data_opt_full[1:10000,8])
deltax_highlight = [synch_error_data_opt_full_active_coupling[k] ? deltax_full[k] : NaN for k in 1:10000]


deltax_alignment_full = abs.(synch_error_data_opt_full[1:10000,5])
deltax_alignment_highlight = [synch_error_data_opt_full_active_coupling[k] ? deltax_alignment_full[k] : NaN for k in 1:10000]


mask = .!isnan.(deltax_alignment_highlight[1:10000]) .& (deltax_alignment_highlight[1:10000] .< 0.5)
# Pad the mask with false to detect edges at the boundaries
padded = [false; mask; false]
# A block starts if current is true and previous was false
starts = findall(i -> padded[i] && !padded[i-1], 2:length(padded))
# A block ends if current is true and next is false
ends = findall(i -> padded[i] && !padded[i+1], 2:length(padded))
datablocks = collect(zip(starts, ends))



#panel c: mechanism

axis_panel_c = 
Axis(
    panel_c[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,22.5 ,1.0e-2, 5e1),
    xticks = [0,5,10,15,20],
    xtickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-2, 1e-1 , 1e0, 1e1],
    ytickformat = values -> [L"10^{-2}",L"10^{-1}",L"10^{0}",L"10^{1}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv"), '\t', Float64)
lines!(
    axis_panel_c,
    synch_error_data_uncoupled[:,1],exp.(synch_error_data_uncoupled[:,8]),
    color = uncoupled_color,
    linewidth = 1.5
)


lines!(
    axis_panel_c,
    t_full, deltax_full,
    color = default_color,
    linewidth = 1.5
)


scatter!(
    axis_panel_c,
    t_full, deltax_highlight,
    color = optimal_coupling_region_timeseries_color,
    markersize = 3.0
)

for (s,e) in datablocks
    vspan!(
        axis_panel_c, 
        t_full[s], 
        t_full[e], 
        color = optimal_coupling_region_timeseries_color,
        alpha = 0.2
    )
end


text!(axis_panel_c, 
    (7, 0.5e1), 
    text = L"\textbf{A}", 
    #space = :relative,
    align = (:center, :middle),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_c,
    (8, 1.2e1),
    (5, 0.0),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    #space = :relative,
)


text!(axis_panel_c, 
    (5, 2e-2), 
    text = L"\textbf{B,C}", 
    #space = :relative,
    align = (:center, :middle),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_c,
    (5.0, 1e-1),
    (-1.5, 7e-1),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    #space = :relative,
)

arrows2d!(axis_panel_c,
    (5.0, 1e-1),
    (1.5, 7e-1),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    #space = :relative,
)




Label(
    panel_c[1, 1, TopLeft()], 
    "c",
    fontsize = 12,
    font = :bold,
    padding = (-8*unit_per_mm, 10*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel d: alignment

#panel l: alignment during coupling

axis_panel_d = 
Axis(
    panel_d[1, 1],
    #ylabel = L"\left|\Delta x\,\right|",
    ylabel = L"\left|\Delta \mathbf{x}\cdot\mathbf{e}_x\,\right|",
    limits = (0,22.5 ,0,1),
    xticks = [0,5,10,15,20],
    xtickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 2.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [0, 0.5, 1.0],
    ytickformat = values -> [L"0", L"0.5", L"1.0"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

lines!(
    axis_panel_d,
    t_full, deltax_alignment_full,
    color = default_color,
    linewidth = 1.5
)

scatter!(
    axis_panel_d,
    t_full, deltax_alignment_highlight,
    color = optimal_coupling_region_timeseries_color,
    markersize = 3.0
)


for (s,e) in datablocks
    vspan!(
        axis_panel_d, 
        t_full[s], 
        t_full[e], 
        color = optimal_coupling_region_timeseries_color,
        alpha = 0.2
    )
end


Label(
    panel_d[1, 1, TopLeft()], 
    "d",
    fontsize = 12,
    font = :bold,
    padding = (-8*unit_per_mm, 10*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)







#panel e: nonlinear simulation

axis_panel_e = 
Axis(
    panel_e[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0, 2.5e5 ,1e-13, 1e2),
    xticks = [0, 1e5, 2e5],
    xtickformat = values -> [L"0", L"10^{5}", L"2\times10^{5}"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-12, 1e-9, 1e-6, 1e-3 , 1e0],
    ytickformat = values -> [L"10^{-12}", L"10^{-9}",L"10^{-6}",L"10^{-3}",L"10^{0}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_uncoupled_nonlinear = DelimitedFiles.readdlm(DrWatson.datadir("roessler","nonlinear__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv"), '\t', Float64)
lines!(
    axis_panel_e,
    synch_error_data_uncoupled_nonlinear[:,1],
    sqrt.( (synch_error_data_uncoupled_nonlinear[:,2] - synch_error_data_uncoupled_nonlinear[:,5]).^2 
            .+ (synch_error_data_uncoupled_nonlinear[:,3] - synch_error_data_uncoupled_nonlinear[:,6]).^2 
            .+ (synch_error_data_uncoupled_nonlinear[:,4] - synch_error_data_uncoupled_nonlinear[:,7]).^2 
           ),
    color = uncoupled_color,
    linewidth = 1.5
)


synch_error_data_opt_full_nonlinear = DelimitedFiles.readdlm(DrWatson.datadir("roessler","nonlinear__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_e,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_opt_full_nonlinear[:,1],
        sqrt.( (synch_error_data_opt_full_nonlinear[:,2] - synch_error_data_opt_full_nonlinear[:,5]).^2 
            .+ (synch_error_data_opt_full_nonlinear[:,3] - synch_error_data_opt_full_nonlinear[:,6]).^2 
            .+ (synch_error_data_opt_full_nonlinear[:,4] - synch_error_data_opt_full_nonlinear[:,7]).^2 
           ),
    color = default_color,
    linewidth = 1.5
)


Label(
    panel_e[1, 1, TopLeft()], 
    "e",
    fontsize = 12,
    font = :bold,
    padding = (-10*unit_per_mm, 10*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)



#panel f: lyap scan for alpha = -1 optimal coupling regions (N = 15)

axis_panel_f = 
Axis(
    panel_f[1, 1],
    ylabel = L"\lambda",
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (-3, 7 ,-0.06, 0.16),
    #xticks = [0,200,400,600,800],
    #xtickformat = values -> [L"0","",L"400","",L"800",""],
    xticks = [-2,0,2,4,6,8,10],
    xtickformat = values -> [L"-2",L"0",L"2", L"4", L"6", L"8"],
    xlabel = L"\text{couping strength}\;\alpha",
    xlabelpadding = -1,
    ylabelpadding = -2.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    #yticks = [1e-18,1e-12,1e-6,1e0,1e6,1e12,1e18,1e24],
    #ytickformat = values -> ["",L"10^{-12}","",L"10^0","",L"10^{12}","",L"10^{24}"],
    yticks = [-0.05,0,0.05,0.10,0.15],
    ytickformat = values -> [L"-0.05", L"0", L"0.05", L"0.10", L"0.15"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)


lines!(
    axis_panel_f,
    [-3,7],[0,0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

lines!(
    axis_panel_f,
    [0,0],[-0.06,0.16],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

lyap_scan_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","lyap_scan__fullycoupled__alpha_1.00_0.00_0.00.tsv"), '\t', Float64)
lines!(
    axis_panel_f,
    lyap_scan_data[:,1],lyap_scan_data[:,4],
    color = default_color,
    alpha = 0.5,
    linewidth = 1.5
)


lyap_scan_opt_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","lyap_scan__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_1.00_0.00_0.00_ncr_15.tsv"), '\t', Float64)
lines!(
    axis_panel_f,
    lyap_scan_opt_data[:,1],lyap_scan_opt_data[:,4],
    color = TUD_darkblue,
    linewidth = 1.5
)

scatter!(
    axis_panel_f,
    [-1],[-0.012013082158072463],
    color = TUD_darkblue,
    markersize = 10
)



Label(
    panel_f[1, 1, TopLeft()], 
    "f",
    fontsize = 12,
    font = :bold,
    padding = (-10*unit_per_mm, 10*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)



#panel g: achieved lyapunov exponent as a function of N or of coupling fraction

axis_panel_g = 
Axis(
    panel_g[1, 1],
    ylabel = L"\lambda",
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0, 0.175 ,-0.06, 0.11),
    #xticks = [0,200,400,600,800],
    #xtickformat = values -> [L"0","",L"400","",L"800",""],
    xticks = [0, 0.05, 0.10, 0.15],
    xtickformat = values -> [L"0",L"0.05",L"0.10", L"0.15"],
    #xlabel = L"\text{couping fraction}\;f_\mathrm{tot}",
    xlabel = L"\text{couping fraction}\;F",
    xlabelpadding = -1,
    ylabelpadding = -2.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    #yticks = [1e-18,1e-12,1e-6,1e0,1e6,1e12,1e18,1e24],
    #ytickformat = values -> ["",L"10^{-12}","",L"10^0","",L"10^{12}","",L"10^{24}"],
    yticks = [-0.05,0,0.05,0.10],
    ytickformat = values -> [L"-0.05", L"0", L"0.05", L"0.10"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)


lines!(
    axis_panel_g,
    [0,0.175],[0,0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

lyap_scan_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","lyap_scan__fullycoupled__alpha_1.00_0.00_0.00.tsv"), '\t', Float64)
lines!(
    axis_panel_g,
    lyap_scan_data[:,1],lyap_scan_data[:,4],
    color = default_color,
    alpha = 0.5,
    linewidth = 1.5
)


lines!(
    axis_panel_g,
    vcat([0.0],optimal_coupling_regions_data[:,9]),vcat([lyap_scan_data[findfirst(==(0.0), lyap_scan_data[:,1]),4]],optimal_coupling_regions_data[:,6]),
    color = TUD_darkblue,
    linewidth = 1.5
)

scatter!(
    axis_panel_g,
    [0.1076589234107659],[-0.012013082158072463],
    color = TUD_darkblue,
    markersize = 10
)



Label(
    panel_g[1, 1, TopLeft()], 
    "g",
    fontsize = 12,
    font = :bold,
    padding = (-10*unit_per_mm, 10*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)





save(plotsdir("FIG2_optimal_coupling_alpha_-1.00_v2.pdf"), fig)
save(plotsdir("FIG2_optimal_coupling_alpha_-1.00_v2.png"), fig)


fig
