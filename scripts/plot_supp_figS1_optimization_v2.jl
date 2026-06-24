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
    figure_padding=3*unit_per_mm,  #1mm padding
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


coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","couplingregions_equidist__dist_0.50__cfrac_0.05.tsv"), '\t', Float64)
optimal_coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00.tsv"), '\t', Float64)


coupling_region_example_1a_index = 482 
coupling_region_example_2a_index = 556

coupling_region_example_1b_index = 482 
coupling_region_example_2b_index = 556

default_color = TUD_darkblue
uncoupled_color = TUD_grey60
optimal_coupling_region_color = TUD_green1
optimal_coupling_region_timeseries_color = RGB(70,178,125)
coupling_region_example_1_color = TUD_orange1
coupling_region_example_2_color = TUD_orange2
coupling_region_example_3_color = TUD_red1

attractor_plot_range = 1:50000
number_of_optimal_coupling_regions = 10



systemparameters_opt_a = SystemParameters(
    StaticArrays.@SVector[5.0, 0.0, 0.0], 
    1, 
    [StaticArrays.@SVector[
        optimal_coupling_regions_data[1,1],
        optimal_coupling_regions_data[1,2],
        optimal_coupling_regions_data[1,3]]
    ], 
    [optimal_coupling_regions_data[1,4]]
)
systemparameters_example_1a = SystemParameters(
    StaticArrays.@SVector[5.0, 0.0, 0.0], 
    1, 
    [StaticArrays.@SVector[
        coupling_regions_data[coupling_region_example_1a_index,1],
        coupling_regions_data[coupling_region_example_1a_index,2],
        coupling_regions_data[coupling_region_example_1a_index,3]]
    ], 
    [coupling_regions_data[coupling_region_example_1a_index,4]]
)
systemparameters_example_2a = SystemParameters(
    StaticArrays.@SVector[5.0, 0.0, 0.0], 
    1, 
    [StaticArrays.@SVector[
        coupling_regions_data[coupling_region_example_2a_index,1],
        coupling_regions_data[coupling_region_example_2a_index,2],
        coupling_regions_data[coupling_region_example_2a_index,3]]
    ], 
    [coupling_regions_data[coupling_region_example_2a_index,4]]
)

simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_opt_a; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__opt_a.tsv", sampling_timestep = 0.01)
simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_example_1a; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_1a.tsv", sampling_timestep = 0.01)
simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_example_2a; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_2a.tsv", sampling_timestep = 0.01)

#lyap_opta = compute_maximum_lyapunov_exponent_linearized(roesslerODE_coupled_linear,"roessler",systemparameters_opt_a)
#lyap_example_1a = compute_maximum_lyapunov_exponent_linearized(roesslerODE_coupled_linear,"roessler",systemparameters_example_1a)
#lyap_example_2a = compute_maximum_lyapunov_exponent_linearized(roesslerODE_coupled_linear,"roessler",systemparameters_example_2a)
#
#print(lyap_opta,"\n")
#print(lyap_example_1a,"\n")
#print(lyap_example_2a,"\n")



systemparameters_opt_b = SystemParameters(
    StaticArrays.@SVector[5.0, 0.0, 0.0], 
    2, 
    [
        StaticArrays.@SVector[
        optimal_coupling_regions_data[1,1],
        optimal_coupling_regions_data[1,2],
        optimal_coupling_regions_data[1,3]],
        StaticArrays.@SVector[
        optimal_coupling_regions_data[2,1],
        optimal_coupling_regions_data[2,2],
        optimal_coupling_regions_data[2,3]]
    ], 
    [
        optimal_coupling_regions_data[1,4],
        optimal_coupling_regions_data[2,4]
    ]
)
systemparameters_example_1b = SystemParameters(
    StaticArrays.@SVector[5.0, 0.0, 0.0], 
    2, 
    [
        StaticArrays.@SVector[
        optimal_coupling_regions_data[1,1],
        optimal_coupling_regions_data[1,2],
        optimal_coupling_regions_data[1,3]],
        StaticArrays.@SVector[
        coupling_regions_data[coupling_region_example_1b_index,1],
        coupling_regions_data[coupling_region_example_1b_index,2],
        coupling_regions_data[coupling_region_example_1b_index,3]]
    ], 
    [
        optimal_coupling_regions_data[1,4],
        coupling_regions_data[coupling_region_example_1b_index,4]
    ]
)
systemparameters_example_2b = SystemParameters(
    StaticArrays.@SVector[5.0, 0.0, 0.0], 
    2, 
    [
        StaticArrays.@SVector[
        optimal_coupling_regions_data[1,1],
        optimal_coupling_regions_data[1,2],
        optimal_coupling_regions_data[1,3]],
        StaticArrays.@SVector[
        coupling_regions_data[coupling_region_example_2b_index,1],
        coupling_regions_data[coupling_region_example_2b_index,2],
        coupling_regions_data[coupling_region_example_2b_index,3]]
    ], 
    [
        optimal_coupling_regions_data[1,4],
        coupling_regions_data[coupling_region_example_2b_index,4]
    ]
)

simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_opt_b; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__opt_b.tsv", sampling_timestep = 0.01)
simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_example_1b; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_1b.tsv", sampling_timestep = 0.01)
simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_example_2b; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_2b.tsv", sampling_timestep = 0.01)

#lyap_optb = compute_maximum_lyapunov_exponent_linearized(roesslerODE_coupled_linear,"roessler",systemparameters_opt_b)
#lyap_example_1b = compute_maximum_lyapunov_exponent_linearized(roesslerODE_coupled_linear,"roessler",systemparameters_example_1b)
#lyap_example_2b = compute_maximum_lyapunov_exponent_linearized(roesslerODE_coupled_linear,"roessler",systemparameters_example_2b)
#
#print(lyap_optb,"\n")
#print(lyap_example_1b,"\n")
#print(lyap_example_2b,"\n")



synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00_ncr_10__trajectory_plot.tsv"), '\t', Float64)
init_state = synch_error_data_opt_full[1,2:8]

systemparameters_uncoupled = SystemParameters(
    StaticArrays.@SVector[0.0, 0.0, 0.0], 
    0, 
    StaticArrays.SVector{3, Float64}[], 
    Float64[]
)

#simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4], init_state[5], init_state[6]], transient_time = 0.0, sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 0.01)
simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4], init_state[5], init_state[6]], transient_time = 0.0, sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 0.05)


#structure of figure

fig = Figure(
            #image size in Makie units converted with 300dpi from 183mm x 90mm (using Makies 0.75 pt_per_unit)
            #size = (ceil(183*25.4*300/0.75), ceil(90*25.4*300/0.75))
        )

#overall figure layout
group_fig = fig[1,1] = GridLayout()

group_abcdefgh = group_fig[1,1] = GridLayout()

group_abcd = group_abcdefgh[1,1] = GridLayout()
panel_a = group_abcd[1,1] = GridLayout()
group_bcd = group_abcd[1,2] = GridLayout()
panel_b = group_bcd[1,1] = GridLayout()
panel_c = group_bcd[2,1] = GridLayout()
panel_d = group_bcd[3,1] = GridLayout()

group_efgh = group_abcdefgh[2,1] = GridLayout()
panel_e = group_efgh[1,1] = GridLayout()
group_fgh = group_efgh[1,2] = GridLayout()
panel_f = group_fgh[1,1] = GridLayout()
panel_g = group_fgh[2,1] = GridLayout()
panel_h = group_fgh[3,1] = GridLayout()

group_ijklm = group_fig[1,2] = GridLayout()
panel_i = group_ijklm[1,1] = GridLayout()
group_jklm = group_ijklm[1,2] = GridLayout()
panel_j = group_jklm[1,1] = GridLayout()
panel_k = group_jklm[2,1] = GridLayout()
panel_l = group_jklm[3,1] = GridLayout()
panel_m = group_jklm[4,1] = GridLayout()



colsize!(group_fig,1,Auto(0.4))
colsize!(group_fig,2,Auto(0.6))

rowsize!(group_abcdefgh,1,Auto(0.5))
rowsize!(group_abcdefgh,2,Auto(0.5))

colsize!(group_abcd,1,Auto(0.66))
colsize!(group_abcd,2,Auto(0.34))

colsize!(group_efgh,1,Auto(0.66))
colsize!(group_efgh,2,Auto(0.34))

colsize!(group_ijklm,1,Auto(0.7))
colsize!(group_ijklm,2,Auto(0.3))



#group_efgh.alignmode = Mixed(top = 0)
#group_ij.alignmode = Mixed(top = 0)


colgap!(group_abcd, 2*unit_per_mm)
colgap!(group_efgh, 2*unit_per_mm)

colgap!(group_ijklm, 2*unit_per_mm)

rowgap!(group_bcd, 1*unit_per_mm)
rowgap!(group_fgh, 1*unit_per_mm)

rowgap!(group_jklm, 2*unit_per_mm)

rowgap!(group_abcdefgh, 3*unit_per_mm)

colgap!(group_fig, 5*unit_per_mm)

#panel_d.alignmode = Mixed(left = 0, right = 0)

#rowsize!(group_abc,1, Auto(0.4))
#rowsize!(group_abc,2, Auto(0.6))

#colsize!(group_ab,1, Auto(0.33))
#colsize!(group_ab,2, Auto(0.66))




#panel_b.alignmode = Mixed(bottom = 0)
#panel_c.alignmode = Mixed(left = 0, right = 0)


#individual panels



#panel a: attractor plot first coupling region

axis_panel_a = 
Axis3(
    panel_a[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.35 * pi,
    elevation = 0.15 * pi,
    limits = (-12,12,-12,8,-3,24),
    xticks = [-10,-5,0,5,10],
    xtickformat = values -> [L"-10","",L"0","",L"10"],
    yticks = [-10,-5,0,5],
    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
    zticks = [0,5,10,15,20],
    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
#    limits = (-10,12,-12,8,-2,25),
#    xticks = [-10,-5,0,5,10],
#    xtickformat = values -> [L"-10","",L"0","",L"10"],
#    yticks = [-10,-5,0,5],
#    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
#    zticks = [0,5,10,15,20,25],
#    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20",L"25"],
    protrusions = (20,0,10,0),
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
    linewidth = 0.25, 
    color = TUD_darkblue,
    alpha = 0.5,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_a, 
    optimal_coupling_regions_data[1,1],
    optimal_coupling_regions_data[1,2],
    optimal_coupling_regions_data[1,3],
    color = coupling_region_example_1_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data[1,1], optimal_coupling_regions_data[1,2], optimal_coupling_regions_data[1,3])]
meshscatter!(axis_panel_a,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data[1,4], 
    color = coupling_region_example_1_color,
    alpha = 0.25,
    transparency=true,
    rasterize = 10
)

#plot other coupling sphere example 1
scatter!(
    axis_panel_a, 
    coupling_regions_data[coupling_region_example_1a_index,1],
    coupling_regions_data[coupling_region_example_1a_index,2],
    coupling_regions_data[coupling_region_example_1a_index,3],
    color = coupling_region_example_2_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(coupling_regions_data[coupling_region_example_1a_index,1], coupling_regions_data[coupling_region_example_1a_index,2], coupling_regions_data[coupling_region_example_1a_index,3])]
meshscatter!(axis_panel_a,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = coupling_regions_data[coupling_region_example_1a_index,4], 
    color = coupling_region_example_2_color,
    alpha = 0.25,
    transparency=true,
    rasterize = 10
)

#plot other coupling sphere example 2
scatter!(
    axis_panel_a, 
    coupling_regions_data[coupling_region_example_2a_index,1],
    coupling_regions_data[coupling_region_example_2a_index,2],
    coupling_regions_data[coupling_region_example_2a_index,3],
    color = coupling_region_example_3_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(coupling_regions_data[coupling_region_example_2a_index,1], coupling_regions_data[coupling_region_example_2a_index,2], coupling_regions_data[coupling_region_example_2a_index,3])]
meshscatter!(axis_panel_a,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = coupling_regions_data[coupling_region_example_2a_index,4], 
    color = coupling_region_example_3_color,
    alpha = 0.25,
    transparency=true,
    rasterize = 10
)

text!(axis_panel_a, 
    (0.085, 1.0), 
    text = L"\alpha = 5.0", 
    space = :relative,
    align = (:left, :top),
    color = :black
    )    

text!(axis_panel_a, 
    (0.975,1.0), 
    text = L"m = 1",
    space = :relative,
    align = (:right, :top),
    color = :black
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



#panel b: example coupling region

axis_panel_b = 
Axis(
    panel_b[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,1e3 ,1e-27, 1e27),
    xticks = [0,200,400,600,800],
    xticklabelsvisible = false,
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-20,1e-10,1e0,1e10,1e20],
    ytickformat = values -> [L"10^{-20}","",L"10^{0}","",L"10^{20}"],
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

synch_error_data_example_1a = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_1a.tsv"), '\t', Float64)
lines!(
    axis_panel_b,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_example_1a[:,1],exp.(synch_error_data_example_1a[:,8]),
    color = coupling_region_example_2_color,
    linewidth = 1.5
)

lines!(
    axis_panel_b,
    [200,900],[1e0,1.86277e18],
    color = coupling_region_example_2_color,
    linewidth = 1
)

text!(
    axis_panel_b,
    (300,1e-15),
    text = L"\lambda \approx 0.060",  #0.060097995
    rotation = 0.08*pi,
    color = coupling_region_example_2_color
)

Label(
    panel_b[1, 1, TopLeft()], 
    "b",
    fontsize = 12,
    font = :bold,
    padding = (11*unit_per_mm, 11*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)



#panel c: example coupling region

axis_panel_c = 
Axis(
    panel_c[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,1e3 ,1e-27, 1e27),
    xticks = [0,200,400,600,800],
    xticklabelsvisible = false,
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-20,1e-10,1e0,1e10,1e20],
    ytickformat = values -> [L"10^{-20}","",L"10^{0}","",L"10^{20}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__uncoupled_trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_c,
    synch_error_data_uncoupled[:,1],exp.(synch_error_data_uncoupled[:,8]),
    color = uncoupled_color,
    linewidth = 1.5
)

synch_error_data_opt_a = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__opt_a.tsv"), '\t', Float64)
lines!(
    axis_panel_c,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_opt_a[:,1],exp.(synch_error_data_opt_a[:,8]),
    color = coupling_region_example_1_color,
    linewidth = 1.5
)

lines!(
    axis_panel_c,
    [200,700],[1e0,1.6290e-17],
    color = coupling_region_example_1_color,
    linewidth = 1
)

text!(
    axis_panel_c,
    (350,1e-3),
    text = L"-0.077",  #-0.077312
    rotation = -0.1*pi,
    color = coupling_region_example_1_color
)

Label(
    panel_c[1, 1, TopLeft()], 
    "c",
    fontsize = 12,
    font = :bold,
    padding = (11*unit_per_mm, 11*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)



#panel d: example coupling region

axis_panel_d = 
Axis(
    panel_d[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,1e3 ,1e-27, 1e27),
    xticks = [0,200,400,600,800],
    xtickformat = values -> [L"0","",L"400","",L"800"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-20,1e-10,1e0,1e10,1e20],
    ytickformat = values -> [L"10^{-20}","",L"10^{0}","",L"10^{20}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__uncoupled_trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_d,
    synch_error_data_uncoupled[:,1],exp.(synch_error_data_uncoupled[:,8]),
    color = uncoupled_color,
    linewidth = 1.5
)

synch_error_data_example_2a = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_2a.tsv"), '\t', Float64)
lines!(
    axis_panel_d,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_example_2a[:,1],exp.(synch_error_data_example_2a[:,8]),
    color = coupling_region_example_3_color,
    linewidth = 1.5
)

lines!(
    axis_panel_d,
    [200,900],[1e-5,3.16e-7],
    color = coupling_region_example_3_color,
    linewidth = 1
)

text!(
    axis_panel_d,
    (250,1e-22),
    text = L"-0.005",  #-0.004935
    rotation = 0.00*pi,
    color = coupling_region_example_3_color
)

Label(
    panel_d[1, 1, TopLeft()], 
    "d",
    fontsize = 12,
    font = :bold,
    padding = (11*unit_per_mm, 11*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)






#panel e: attractor plot second coupling region

axis_panel_e = 
Axis3(
    panel_e[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.35 * pi,
    elevation = 0.15 * pi,
    limits = (-12,12,-12,8,-3,24),
    xticks = [-10,-5,0,5,10],
    xtickformat = values -> [L"-10","",L"0","",L"10"],
    yticks = [-10,-5,0,5],
    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
    zticks = [0,5,10,15,20],
    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
#    limits = (-10,12,-12,8,-2,25),
#    xticks = [-10,-5,0,5,10],
#    xtickformat = values -> [L"-10","",L"0","",L"10"],
#    yticks = [-10,-5,0,5],
#    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
#    zticks = [0,5,10,15,20,25],
#    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20",L"25"],
    protrusions = (20,0,10,0),
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


lines!(axis_panel_e, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.25, 
    color = TUD_darkblue,
    alpha = 0.5,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_e, 
    optimal_coupling_regions_data[1,1],
    optimal_coupling_regions_data[1,2],
    optimal_coupling_regions_data[1,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data[1,1], optimal_coupling_regions_data[1,2], optimal_coupling_regions_data[1,3])]
meshscatter!(axis_panel_e,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data[1,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency=true,
    rasterize = 10
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_e, 
    optimal_coupling_regions_data[2,1],
    optimal_coupling_regions_data[2,2],
    optimal_coupling_regions_data[2,3],
    color = coupling_region_example_1_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data[2,1], optimal_coupling_regions_data[2,2], optimal_coupling_regions_data[2,3])]
meshscatter!(axis_panel_e,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data[2,4], 
    color = coupling_region_example_1_color,
    alpha = 0.25,
    transparency=true,
    rasterize = 10
)

#plot other coupling sphere example 1
scatter!(
    axis_panel_e, 
    coupling_regions_data[coupling_region_example_1b_index,1],
    coupling_regions_data[coupling_region_example_1b_index,2],
    coupling_regions_data[coupling_region_example_1b_index,3],
    color = coupling_region_example_2_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(coupling_regions_data[coupling_region_example_1b_index,1], coupling_regions_data[coupling_region_example_1b_index,2], coupling_regions_data[coupling_region_example_1b_index,3])]
meshscatter!(axis_panel_e,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = coupling_regions_data[coupling_region_example_1b_index,4], 
    color = coupling_region_example_2_color,
    alpha = 0.25,
    transparency=true,
    rasterize = 10
)

#plot other coupling sphere example 2
scatter!(
    axis_panel_e, 
    coupling_regions_data[coupling_region_example_2b_index,1],
    coupling_regions_data[coupling_region_example_2b_index,2],
    coupling_regions_data[coupling_region_example_2b_index,3],
    color = coupling_region_example_3_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(coupling_regions_data[coupling_region_example_2b_index,1], coupling_regions_data[coupling_region_example_2b_index,2], coupling_regions_data[coupling_region_example_2b_index,3])]
meshscatter!(axis_panel_e,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = coupling_regions_data[coupling_region_example_2b_index,4], 
    color = coupling_region_example_3_color,
    alpha = 0.25,
    transparency=true,
    rasterize = 10
)

text!(axis_panel_e, 
    (0.085, 1.0), 
    text = L"\alpha = 5.0", 
    space = :relative,
    align = (:left, :top),
    color = :black
    )    

text!(axis_panel_e, 
    (0.975,1.0), 
    text = L"m = 2", 
    space = :relative,
    align = (:right, :top),
    color = :black
    )  

Label(
    panel_e[1, 1, TopLeft()], 
    "e",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)



#panel f: example coupling region

axis_panel_f = 
Axis(
    panel_f[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,1e3 ,1e-27, 1e6),
    xticks = [0,200,400,600,800],
    xticklabelsvisible = false,
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-20,1e-10,1e0,1e10,1e20],
    ytickformat = values -> [L"10^{-20}","",L"10^{0}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_opt_a = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__opt_a.tsv"), '\t', Float64)
lines!(
    axis_panel_f,
    synch_error_data_opt_a[:,1],exp.(synch_error_data_opt_a[:,8]),
    color = TUD_grey20,
    linewidth = 1.5
)

synch_error_data_example_1b = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_1b.tsv"), '\t', Float64)
lines!(
    axis_panel_f,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_example_1b[:,1],exp.(synch_error_data_example_1b[:,8]),
    color = coupling_region_example_2_color,
    linewidth = 1.5
)

lines!(
    axis_panel_f,
    [250,850],[1e0,9.66085e-14],
    color = coupling_region_example_2_color,
    linewidth = 1
)

text!(
    axis_panel_f,
    (400,2e-5),
    text = L"-0.050",  #-0.04994685
    rotation = -0.09*pi,
    color = coupling_region_example_2_color
)

Label(
    panel_f[1, 1, TopLeft()], 
    "f",
    fontsize = 12,
    font = :bold,
    padding = (11*unit_per_mm, 11*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)



#panel g: example coupling region

axis_panel_g = 
Axis(
    panel_g[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,1e3 ,1e-27, 1e6),
    xticks = [0,200,400,600,800],
    xticklabelsvisible = false,
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-20,1e-10,1e0,1e10,1e20],
    ytickformat = values -> [L"10^{-20}","",L"10^{0}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_opt_a = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__opt_a.tsv"), '\t', Float64)
lines!(
    axis_panel_g,
    synch_error_data_opt_a[:,1],exp.(synch_error_data_opt_a[:,8]),
    color = TUD_grey20,
    linewidth = 1.5
)

synch_error_data_opt_b = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__opt_b.tsv"), '\t', Float64)
lines!(
    axis_panel_g,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_opt_b[:,1],exp.(synch_error_data_opt_b[:,8]),
    color = coupling_region_example_1_color,
    linewidth = 1.5
)

lines!(
    axis_panel_g,
    [85,360],[1e0,2.8092632e-24],
    color = coupling_region_example_1_color,
    linewidth = 1
)

text!(
    axis_panel_g,
    (150,1e0),
    text = L"-0.197",  #-0.1971968
    rotation = -0.33*pi,
    color = coupling_region_example_1_color
)

Label(
    panel_g[1, 1, TopLeft()], 
    "g",
    fontsize = 12,
    font = :bold,
    padding = (11*unit_per_mm, 11*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)



#panel h: example coupling region

axis_panel_h = 
Axis(
    panel_h[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,1e3 ,1e-27, 1e6),
    xticks = [0,200,400,600,800],
    xtickformat = values -> [L"0","",L"400","",L"800"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-20,1e-10,1e0],
    ytickformat = values -> [L"10^{-20}","",L"10^{0}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_opt_a = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__opt_a.tsv"), '\t', Float64)
lines!(
    axis_panel_h,
    synch_error_data_opt_a[:,1],exp.(synch_error_data_opt_a[:,8]),
    color = TUD_grey20,
    linewidth = 1.5
)

synch_error_data_example_2b = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__example_2b.tsv"), '\t', Float64)
lines!(
    axis_panel_h,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_example_2b[:,1],exp.(synch_error_data_example_2b[:,8]),
    color = coupling_region_example_3_color,
    linewidth = 1.5
)

lines!(
    axis_panel_h,
    [200,800],[1e-5,1.14951e-13],
    color = coupling_region_example_3_color,
    linewidth = 1
)

text!(
    axis_panel_h,
    (200,1e-15),
    text = L"-0.030",  #-0.0304689095
    rotation = -0.075*pi,
    color = coupling_region_example_3_color
)

Label(
    panel_h[1, 1, TopLeft()], 
    "h",
    fontsize = 12,
    font = :bold,
    padding = (11*unit_per_mm, 11*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)





#panel i: attractor plot optimal coupling regions

axis_panel_i = 
Axis3(
    panel_i[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.35 * pi,
    elevation = 0.15 * pi,
    limits = (-12,12,-12,8,-3,24),
    xticks = [-10,-5,0,5,10],
    xtickformat = values -> [L"-10","",L"0","",L"10"],
    yticks = [-10,-5,0,5],
    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
    zticks = [0,5,10,15,20],
    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
#    limits = (-10,12,-12,8,-2,25),
#    xticks = [-10,-5,0,5,10],
#    xtickformat = values -> [L"-10","",L"0","",L"10"],
 #   yticks = [-10,-5,0,5],
 #   ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
 #   zticks = [0,5,10,15,20,25],
 #   ztickformat = values -> [L"0",L"5",L"10",L"15",L"20",L"25"],
    protrusions = (20,0,10,0),
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


lines!(axis_panel_i, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00_ncr_10__trajectory_plot.tsv"), '\t', Float64)
arrow_plot_range = 1:150:200000
arrows3d!(
    axis_panel_i,
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
    axis_panel_i, 
    optimal_coupling_regions_data[1:10,1],
    optimal_coupling_regions_data[1:10,2],
    optimal_coupling_regions_data[1:10,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency = true,
    rasterize = 10
)
spherecenter_positive = [Point(optimal_coupling_regions_data[k,1], optimal_coupling_regions_data[k,2], optimal_coupling_regions_data[k,3]) for k in 1:10]
meshscatter!(axis_panel_i,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data[1:10,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)



text!(axis_panel_i, 
    (0.085, 0.96), 
    text = L"\alpha = 5.0", 
    space = :relative, 
    align = (:left, :top), 
    color = :black, 
    fontsize = 14
    )    

text!(axis_panel_i, 
    (0.975, 0.96), 
    text = L"m = 10", 
    space = :relative, 
    align = (:right, :top), 
    color = :black,
    fontsize = 14
    )  



text!(axis_panel_i, 
    (0.455, 0.18), 
    text = L"\textbf{A}", 
    space = :relative,
    align = (:right, :top),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_i,
    (0.465, 0.165),
    (0.1, 0.05),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    space = :relative,
)

text!(axis_panel_i, 
    (0.445, 0.615), 
    text = L"\textbf{B}", 
    space = :relative,
    align = (:right, :top),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_i,
    (0.46, 0.565),
    (0.05, -0.075),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    space = :relative,
)


text!(axis_panel_i, 
    (0.19, 0.49), 
    text = L"\textbf{C}", 
    space = :relative,
    align = (:right, :top),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_i,
    (0.195, 0.445),
    (0.05, -0.05),
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
    panel_i[1, 1, TopLeft()], 
    "i",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)



#panel j: optimal synch error

axis_panel_j = 
Axis(
    panel_j[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,1e2 ,1e-25, 1e5),
    xticks = [0,20,40,60,80],
    xtickformat = values -> [L"0","",L"40","",L"80"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-20,1e-10,1e0],
    ytickformat = values -> [L"10^{-20}",L"10^{-10}",L"10^0"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__uncoupled_trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_j,
    synch_error_data_uncoupled[:,1],exp.(synch_error_data_uncoupled[:,8]),
    color = uncoupled_color,
    linewidth = 1.5
)

synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00_ncr_10__trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_j,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_opt_full[:,1],exp.(synch_error_data_opt_full[:,8]),
    color = TUD_darkblue,
    linewidth = 1.5
)

lines!(
    axis_panel_j,
    [25,90],[2.0813e-2,2.8620044e-24],
    color = TUD_darkblue,
    linewidth = 1
)

text!(
    axis_panel_j,
    (45,1e-8),
    text = L"\lambda \approx -0.774",  #-0.774436
    rotation = -0.17*pi,
    color = TUD_darkblue
)

lines!(
    axis_panel_j,
    [0,1000],[1e0,1e0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

Label(
    panel_j[1, 1, TopLeft()], 
    "j",
    fontsize = 12,
    font = :bold,
    padding = (12*unit_per_mm, 12*unit_per_mm, -3*unit_per_mm, -3*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)




synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00_ncr_10__trajectory_plot.tsv"), '\t', Float64)
synch_error_data_opt_full_active_coupling = [ 
    is_in_coupling_region(synch_error_data_opt_full[k,2:4], "couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00.tsv", number_of_optimal_coupling_regions) 
    for k in 1:10000
]


t_full = synch_error_data_opt_full[1:10000,1]
deltax_full = exp.(synch_error_data_opt_full[1:10000,8])
deltax_highlight = [synch_error_data_opt_full_active_coupling[k] ? deltax_full[k] : NaN for k in 1:10000]


deltax_alignment_full = abs.(synch_error_data_opt_full[1:10000,5])
deltax_alignment_highlight = [synch_error_data_opt_full_active_coupling[k] ? deltax_alignment_full[k] : NaN for k in 1:10000]


mask = .!isnan.(deltax_alignment_highlight[1:10000]) .& (deltax_alignment_highlight[1:10000] .> 0.7)
# Pad the mask with false to detect edges at the boundaries
padded = [false; mask; false]
# A block starts if current is true and previous was false
starts = findall(i -> padded[i] && !padded[i-1], 2:length(padded))
# A block ends if current is true and next is false
ends = findall(i -> padded[i] && !padded[i+1], 2:length(padded))
datablocks = collect(zip(starts, ends))



#panel k: trajectory and mechanism

axis_panel_k = 
Axis(
    panel_k[1, 1],
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0,22.5 ,5e-8, 5e0),
    xticks = [0,5,10,15,20],
    xtickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
    xticklabelsvisible = false,
    xlabelpadding = -1,
    ylabelpadding = 3.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0],
    ytickformat = values -> [L"10^{-6}","","",L"10^{-3}","","",L"10^{0}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

lines!(
    axis_panel_k,
    [0,22.5],[1e0,1e0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)


synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_5.00_0.00_0.00__uncoupled_comparison.tsv"), '\t', Float64)
lines!(
    axis_panel_k,
    synch_error_data_uncoupled[:,1],exp.(synch_error_data_uncoupled[:,8]),
    color = uncoupled_color,
    linewidth = 1.5
)


lines!(
    axis_panel_k,
    t_full, deltax_full,
    color = default_color,
    linewidth = 1.5
)


scatter!(
    axis_panel_k,
    t_full, deltax_highlight,
    color = optimal_coupling_region_timeseries_color,
    markersize = 3.0
)

for (s,e) in datablocks
    vspan!(
        axis_panel_k, 
        t_full[s], 
        t_full[e], 
        color = optimal_coupling_region_timeseries_color,
        alpha = 0.2
    )
end

text!(axis_panel_k, 
    (19, 1e-3), 
    text = L"\textbf{A}", 
    #space = :relative,
    align = (:center, :middle),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_k,
    (19, 1e-3),
    (-2.5, -9e-4),
    shaftcolor = optimal_coupling_region_color, 
    tipcolor = optimal_coupling_region_color,
    shaftwidth = 2.0,
    tipwidth = 10.0,
    taillength = 0.0,
    tailwidth = 0.0,
    lengthscale = 1,
    #space = :relative,
)


text!(axis_panel_k, 
    (13, 0.4e-2), 
    text = L"\textbf{B,C}", 
    #space = :relative,
    align = (:center, :middle),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_k,
    (13.5, 0.5e-2),
    (0, -0.45e-2),
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
    panel_k[1, 1, TopLeft()], 
    "k",
    fontsize = 12,
    font = :bold,
    padding = (12*unit_per_mm, 12*unit_per_mm, -3*unit_per_mm, -3*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)





#panel l: alignment during coupling

axis_panel_l = 
Axis(
    panel_l[1, 1],
    ylabel = L"\left|\Delta \mathbf{x}\cdot\mathbf{e}_x\,\right|",
    limits = (0,22.5 ,0,1),
    xticks = [0,5,10,15,20],
    xtickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 10.5,
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
    axis_panel_l,
    t_full, deltax_alignment_full,
    color = default_color,
    linewidth = 1.5
)

scatter!(
    axis_panel_l,
    t_full, deltax_alignment_highlight,
    color = optimal_coupling_region_timeseries_color,
    markersize = 3.0
)


for (s,e) in datablocks
    vspan!(
        axis_panel_l, 
        t_full[s], 
        t_full[e], 
        color = optimal_coupling_region_timeseries_color,
        alpha = 0.2
    )
end


Label(
    panel_l[1, 1, TopLeft()], 
    "l",
    fontsize = 12,
    font = :bold,
    padding = (12*unit_per_mm, 12*unit_per_mm, -3*unit_per_mm, -3*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)






#panel m: lyap scan for alpha = 5 optimal coupling regions

axis_panel_m = 
Axis(
    panel_m[1, 1],
    ylabel = L"\lambda",
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0, 7 ,-1.1, 0.1),
    #xticks = [0,200,400,600,800],
    #xtickformat = values -> [L"0","",L"400","",L"800",""],
    xticks = [0,2,4,6,8,10],
    xtickformat = values -> [L"0",L"2", L"4", L"6", L"8"],
    xlabel = L"\text{couping strength}\;\alpha",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    #yticks = [1e-18,1e-12,1e-6,1e0,1e6,1e12,1e18,1e24],
    #ytickformat = values -> ["",L"10^{-12}","",L"10^0","",L"10^{12}","",L"10^{24}"],
    yticks = [-1.0, -0.5, 0 ],
    ytickformat = values -> [L"-1.0", L"-0.5", L"0"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)


lines!(
    axis_panel_m,
    [0,10],[0,0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

lyap_scan_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","lyap_scan__fullycoupled__alpha_1.00_0.00_0.00.tsv"), '\t', Float64)
lines!(
    axis_panel_m,
    lyap_scan_data[:,1],lyap_scan_data[:,4],
    color = default_color,
    alpha = 0.5,
    linewidth = 1.5
)


lyap_scan_opt_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","lyap_scan__couplingregions_equidist__dist_0.50__cfrac_0.05__optimal_coupling__alpha_5.00_0.00_0.00_ncr_10.tsv"), '\t', Float64)
lines!(
    axis_panel_m,
    lyap_scan_opt_data[:,1],lyap_scan_opt_data[:,4],
    color = TUD_darkblue,
    linewidth = 1.5
)

scatter!(
    axis_panel_m,
    [5],[-0.7744360700744921],
    color = TUD_darkblue,
    markersize = 10
)



Label(
    panel_m[1, 1, TopLeft()], 
    "m",
    fontsize = 12,
    font = :bold,
    padding = (10*unit_per_mm, 10*unit_per_mm, -3*unit_per_mm, -3*unit_per_mm),
    halign = :right,
    valign = :bottom,
    tellwidth = false,
    tellheight = false
)



save(plotsdir("FIGS1_optimal_coupling_alpha_5.00_v2.pdf"), fig)
save(plotsdir("FIGS1_optimal_coupling_alpha_5.00_v2.png"), fig)


fig