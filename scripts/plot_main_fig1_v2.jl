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
    size=(663,275),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
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
coupling_region_positive_index = 674
coupling_region_negative_index = 855

default_color = TUD_darkblue
uncoupled_color = TUD_grey60
coupling_region_positive_color = TUD_green1
coupling_region_negative_color = TUD_orange1

attractor_plot_range = 1:50000


systemparameters_positive = SystemParameters(
    StaticArrays.@SVector[1.0, 0.0, 0.0], 
    1, 
    [StaticArrays.@SVector[
        coupling_regions_data[coupling_region_positive_index,1],
        coupling_regions_data[coupling_region_positive_index,2],
        coupling_regions_data[coupling_region_positive_index,3]]
    ], 
    [coupling_regions_data[coupling_region_positive_index,4]]
)
systemparameters_negative = SystemParameters(
    StaticArrays.@SVector[1.0, 0.0, 0.0], 
    1, 
    [StaticArrays.@SVector[
        coupling_regions_data[coupling_region_negative_index,1],
        coupling_regions_data[coupling_region_negative_index,2],
        coupling_regions_data[coupling_region_negative_index,3]]
    ], 
    [coupling_regions_data[coupling_region_negative_index,4]]
)

simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_positive; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00__positive_example.tsv", sampling_timestep = 0.01)
simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_negative; sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00__negative_example.tsv", sampling_timestep = 0.01)


#structure of figure

fig = Figure(
            #image size in Makie units converted with 300dpi from 183mm x 90mm (using Makies 0.75 pt_per_unit)
            #size = (ceil(183*25.4*300/0.75), ceil(90*25.4*300/0.75))
        )

#overall figure layout
group_fig = fig[1,1] = GridLayout()
group_abc = group_fig[1,1] = GridLayout()
group_ab = group_abc[1,1] = GridLayout()
panel_a = group_ab[1,1] = GridLayout()
panel_b = group_ab[1,2] = GridLayout()
panel_c = group_abc[2,1] = GridLayout()
group_def = group_fig[1,2] = GridLayout()
panel_d = group_def[1,1] = GridLayout()
group_ef = group_def[1,2] = GridLayout()
panel_e = group_ef[1,1] = GridLayout()
panel_f = group_ef[2,1] = GridLayout()





#colgap!(group_fig,3*unit_per_mm)

#rowgap!(group_efgh, 2*unit_per_mm)

#rowgap!(group_ij, 2*unit_per_mm)

#group_efgh.alignmode = Mixed(top = 0)
#group_ij.alignmode = Mixed(top = 0)




colgap!(group_fig, 1*unit_per_mm)
#rowgap!(group_abcd, 1*unit_per_mm)

rowgap!(group_abc, 3*unit_per_mm)
colgap!(group_ab, -1*unit_per_mm)

colgap!(group_def, 4*unit_per_mm)

rowgap!(group_ef, 3*unit_per_mm)

#panel_d.alignmode = Mixed(left = 0, right = 0)


rowsize!(group_abc,1, Auto(0.4))
rowsize!(group_abc,2, Auto(0.6))

colsize!(group_ab,1, Auto(0.33))
colsize!(group_ab,2, Auto(0.66))


colsize!(group_fig,1,Auto(0.3))
colsize!(group_fig,2,Auto(0.7))
colsize!(group_def,1,Auto(0.45))
colsize!(group_def,2,Auto(0.25))

rowsize!(group_ef,1, Auto(0.575))
rowsize!(group_ef,2, Auto(0.425))

panel_b.alignmode = Mixed(bottom = 0)
panel_c.alignmode = Mixed(left = 0, right = 0)


#individual panels

#panel a: master-slave illustration

axis_panel_a =
Axis(
    panel_a[1,1],
    limits = (-1.75,1.75,-3,3),
    aspect = 3.5/6
)
hidedecorations!(axis_panel_a)
hidespines!(axis_panel_a)

top_y = 2.0      
bottom_y = -2.0  
r = 0.75    
arrow_y_start = top_y - r - 0.1
arrow_y_end = bottom_y + r + 0.1
arrow_len = arrow_y_end - arrow_y_start


poly!(axis_panel_a, Circle(Point2f(0, top_y), r), color = :white, strokecolor = :black, strokewidth = 2)
poly!(axis_panel_a, Circle(Point2f(0, bottom_y), r), color = :white, strokecolor = :black, strokewidth = 2)


text!(axis_panel_a, 0, top_y; text = L"\mathbf{x}_1", align = (:center, :center), fontsize = 10)
text!(axis_panel_a, 0, bottom_y; text = L"\mathbf{x}_2", align = (:center, :center), fontsize = 10)

arrows2d!(axis_panel_a, [0], [arrow_y_start], [0], [arrow_len], 
        color = :black, shaftwidth = 2, tipwidth = 6)

text!(axis_panel_a, -0.75, (top_y + bottom_y)/2+0.2; 
      text = L"\alpha", align = (:left, :center), fontsize = 10)

Label(
    panel_a[1, 1, TopLeft()], 
    "a",
    fontsize = 12,
    #color = TUD_grey40,
    font = :bold,
    padding = (-2*unit_per_mm, -3*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel b: absolute synchronization difference trajectory data (for alpha = 0, alpha = 1, and alpha = 5)

axis_panel_b = 
Axis(
    panel_b[1, 1],
    xlabel = L"\text{time}\;t",
    ylabel = L"\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    limits = (0,1e2 ,1e-18, 1e8),
    xticks = [0,50,100],
    xtickformat = values -> [L"0", L"50", L"100\text{    }"],
    #xticklabelalign = (:right,:top),
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1,
    yticks = [1e-18,1e-12,1e-6,1e0,1e6],
    ytickformat = values -> [L"10^{-18}",L"10^{-12}",L"10^{-6}",L"10^0",L"10^6"],
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

lines!(
    axis_panel_b,
    [10,90],[1e2,1*300*1e2],
    color = uncoupled_color,
    linewidth = 1
)

text!(
    axis_panel_b,
    (15,0.5e3),
    text = L"\lambda_0 \approx 0.071",  #0.0711981321473530
    rotation = 0.025*pi,
    color = uncoupled_color
)


lines!(
    axis_panel_b,
    [0,100],[1e0,1e0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)


synch_error_data_fullycoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00.tsv"), '\t', Float64)
lines!(
    axis_panel_b,
    synch_error_data_fullycoupled[:,1],exp.(synch_error_data_fullycoupled[:,8]),
    color = default_color,
    linewidth = 1.5
)

lines!(
    axis_panel_b,
    [10,80],[1e-5,1.463e-16],
    color = default_color,
    linewidth = 1
)

text!(
    axis_panel_b,
    (15,4e-11),
    text = L"\lambda \approx -0.356",  #-0.35640052505204933
    rotation = 1.87*pi,
    color = default_color
)




#
# add proportional lines and approximate lyap exponents
#

Label(
    panel_b[1, 1, TopLeft()], 
    "b",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 8*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel c: lyapunov exponent as function of coupling strength

axis_panel_c = 
Axis(
    panel_c[1, 1],
    xlabel = L"\text{coupling strength}\;\alpha",
    ylabel = L"\text{max. Lyap. exp.}\;\;\lambda",
    limits = (0,5,-0.6,0.1),
    xticks = [0,1,2,3,4,5],
    xtickformat = values -> [L"0",L"1",L"2",L"3",L"4",L"5"],
    xlabelpadding = -1,
    ylabelpadding = 2.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    yticks = [-0.6,-0.4,-0.2,0],
    ytickformat = values -> [L"-0.6",L"-0.4",L"-0.2",L"0"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

lines!(
    axis_panel_c,
    [0,5],[0,0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)


#
# add proportional lines and approximate lyap exponents
#
scatter!(
    axis_panel_c,
    [0],[0.07119813214735302],
    color = uncoupled_color,
    markersize = 10
)


scatter!(
    axis_panel_c,
    [1],[-0.35640052505204933],
    color = default_color,
    markersize = 10
)



lyap_scan_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","lyap_scan__fullycoupled__alpha_1.00_0.00_0.00.tsv"), '\t', Float64)
lines!(
    axis_panel_c,
    lyap_scan_data[:,1],lyap_scan_data[:,4],
    color = default_color,
    linewidth = 1.5
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




#panel d: attractor plot 

axis_panel_d = 
Axis3(
    panel_d[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.35 * pi,
    elevation = 0.15 * pi,
    limits = (-10,12,-12,8,-2,25),
    xticks = [-10,-5,0,5,10],
    xtickformat = values -> [L"-10",L"-5",L"0",L"5",L"10"],
    yticks = [-10,-5,0,5],
    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
    zticks = [0,5,10,15,20,25],
    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20",L"25"],
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


lines!(axis_panel_d, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.5,
    transparency = true
)

couplingregion_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","couplingregions_equidist__dist_0.50__cfrac_0.01__lyap_effect__alpha_1.00_0.00_0.00.tsv"), '\t', Float64)
sc = scatter!(
    axis_panel_d, 
    couplingregion_data[:,1],
    couplingregion_data[:,2],
    couplingregion_data[:,3],
    color = couplingregion_data[:,9], 
    colormap = :berlin,
    colorscale = Makie.Symlog10(0.001),
    colorrange = (-0.25, 0.25),
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     
    alpha = 1.0,
    transparency=true
)

# 2. Add the Colorbar in the next column of the grid layout
Colorbar(
    panel_d[1, 2], 
    sc, 
    label = L"\lambda - \lambda_0", # Replace with your actual variable label
    width = 8,
    height = Relative(0.725),
    halign = :center,
    tellheight = false,
    labelpadding = -15,
    ticklabelpad = 2,
    ticksize = 4,
    tickalign = 1,
    scale = Makie.Symlog10(0.001), # Match the scale of your scatter plot
    ticks = ([-1e-1, -1e-2, -1e-3, 0, 1e-3, 1e-2, 1e-1],  # The data values where ticks should appear
             [L"-10^{-1}", L"-10^{-2}", L"-10^{-3}", L"0", L"10^{-3}", L"10^{-2}", L"10^{-1}"]), # The LaTeX labels
)

# Optional: Adjust the spacing between the axis and colorbar
colgap!(panel_d, 1*unit_per_mm)

text!(axis_panel_d, 
    (0.085, 1.00), 
    text = L"\alpha = 1.0", 
    space = :relative,
    align = (:left, :top),
    color = :black,
    fontsize = 14
    )    

#text!(axis_panel_d, 
#    (0.975, 0.925), 
#    text = L"N = 1", 
#    space = :relative,
#    align = (:right, :top),
#    color = :black
#    ) 

Label(
    panel_d[1, 1, TopLeft()], 
    "d",
    fontsize = 12,
    font = :bold,
    padding = (4*unit_per_mm, 5*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel e: example of coupling on synch error evolution (positive coupling -> worse)

axis_panel_e = 
Axis3(
    panel_e[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.35 * pi,
    elevation = 0.15 * pi,
    limits = (-10,12,-12,8,-2,25),
    xticks = [-10,-5,0,5,10],
    xtickformat = values -> [L"-10","",L"0","",L"10"],
    yticks = [-10,-5,0,5],
    ytickformat = values -> [L"-10",L"-5",L"0",L"5"],
    zticks = [0,5,10,15,20,25],
    ztickformat = values -> [L"0",L"5",L"10",L"15",L"20",L"25"],
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

    #plot positive coupling sphere example
    scatter!(
        axis_panel_e, 
        coupling_regions_data[coupling_region_positive_index,1],
        coupling_regions_data[coupling_region_positive_index,2],
        coupling_regions_data[coupling_region_positive_index,3],
        color = coupling_region_positive_color, 
        markersize = 5,      # Set a clear marker size
        strokewidth = 0.0,     # Remove stroke for cleaner appearance
        alpha = 1.0,
        transparency=true
    )
    spherecenter_positive = [Point(coupling_regions_data[coupling_region_positive_index,1], coupling_regions_data[coupling_region_positive_index,2], coupling_regions_data[coupling_region_positive_index,3])]
    meshscatter!(axis_panel_e,
        spherecenter_positive,
        marker = Sphere(Point(0.0,0.0,0.0),1),
        markersize = coupling_regions_data[coupling_region_positive_index,4], 
        color = coupling_region_positive_color,
        alpha = 0.25,
        transparency=true,
        rasterize = 10
    )


    #plot negative coupling sphere example
    scatter!(
        axis_panel_e, 
        coupling_regions_data[coupling_region_negative_index,1],
        coupling_regions_data[coupling_region_negative_index,2],
        coupling_regions_data[coupling_region_negative_index,3],
        color = coupling_region_negative_color, 
        markersize = 5,      # Set a clear marker size
        strokewidth = 0.0,     # Remove stroke for cleaner appearance
        alpha = 1.0,
        transparency=true
    )
    spherecenter_positive = [Point(coupling_regions_data[coupling_region_negative_index,1], coupling_regions_data[coupling_region_negative_index,2], coupling_regions_data[coupling_region_negative_index,3])]
    meshscatter!(axis_panel_e,
        spherecenter_positive,
        marker = Sphere(Point(0.0,0.0,0.0),1),
        markersize = coupling_regions_data[coupling_region_negative_index,4], 
        color = coupling_region_negative_color,
        alpha = 0.25,
        transparency=true,
        rasterize = 10
    )

Label(
    panel_e[1, 1, TopLeft()], 
    "e",
    fontsize = 12,
    font = :bold,
    padding = (-7*unit_per_mm, 10*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel f: example of coupling on synch error evolution (positive coupling -> better)

axis_panel_f = 
Axis(
    panel_f[1, 1],
    xlabel = L"\text{time}\;t",
    ylabel = L"\text{synch. error}\;\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    limits = (0,1e3 ,1e-1, 1e25),
    xticks = [0,500,1000],
    xtickformat = values -> [L"0", L"500", L"1000\text{     }"],
    #xticklabelalign = (:right,:top),
    xlabelpadding = -1,
    ylabelpadding = 2.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    yticks = [1e0,1e6,1e12,1e18,1e24],
    ytickformat = values -> [L"10^0",L"10^{6}",L"10^{12}",L"10^{18}",L"10^{24}"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)


coupling_regions_data = DelimitedFiles.readdlm(DrWatson.datadir("roessler","couplingregions_equidist__dist_0.50__cfrac_0.01.tsv"), '\t', Float64)
for k in 1:size(coupling_regions_data,1)
#for k in 1:50

    systemparameters_temp = SystemParameters(
        StaticArrays.@SVector[1.0, 0.0, 0.0], 
        1, 
        [StaticArrays.@SVector[
            coupling_regions_data[k,1],
            coupling_regions_data[k,2],
            coupling_regions_data[k,3]]
        ], 
        [coupling_regions_data[k,4]]
    )

    simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_temp; sampling = true, output_filename = "temp_sim_data.tsv", sampling_timestep = 4.0, simulation_time = 1010.0)

    synch_error_data_temp = DelimitedFiles.readdlm(DrWatson.datadir("roessler","temp_sim_data.tsv"), '\t', Float64)
    lines!(
        axis_panel_f,
        synch_error_data_temp[:,1],exp.(synch_error_data_temp[:,8]),
        color = uncoupled_color,
        alpha = 0.01,
        #alpha = 0.1,
        linewidth = 0.5
    )
end


synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__uncoupled_trajectory_plot.tsv"), '\t', Float64)
lines!(
    axis_panel_f,
    synch_error_data_uncoupled[:,1],exp.(synch_error_data_uncoupled[:,8]),
    color = uncoupled_color,
    linewidth = 1.5
)

synch_error_data_coupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00__negative_example.tsv"), '\t', Float64)
lines!(
    axis_panel_f,
    #synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8].-synch_error_data_uncoupled[:,8]),
    synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8]),
    color = coupling_region_negative_color,
    linewidth = 1.5
)

lines!(
    axis_panel_f,
    [100,600],[40*5e5,40*22.65e21],
    color = coupling_region_negative_color,
    linewidth = 1
)

text!(
    axis_panel_f,
    (200,2e11),
    text = L"\lambda \approx 0.077",  #0.077....
    rotation = 0.225*pi,
    color = coupling_region_negative_color
)

synch_error_data_coupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_1.00_0.00_0.00__positive_example.tsv"), '\t', Float64)
lines!(
    axis_panel_f,
    synch_error_data_coupled[:,1],exp.(synch_error_data_coupled[:,8]),
    color = coupling_region_positive_color,
    linewidth = 1.5
)

lines!(
    axis_panel_f,
    [200,900],[0.2e2,2.73e15*0.2],
    color = coupling_region_positive_color,
    linewidth = 1
)

text!(
    axis_panel_f,
    (325,10e-1),
    text = L"\lambda \approx 0.044",  #0.04419579781401635
    rotation = 0.175*pi,
    color = coupling_region_positive_color
)

Label(
    panel_f[1, 1, TopLeft()], 
    "f",
    fontsize = 12,
    font = :bold,
    padding = (-7*unit_per_mm, 10*unit_per_mm, -1*unit_per_mm, -3*unit_per_mm),
    halign = :right,
    valign = :bottom
)


save(plotsdir("FIG1_state_dependent_coupling_v2.pdf"), fig)
save(plotsdir("FIG1_state_dependent_coupling_v2.png"), fig)


fig
