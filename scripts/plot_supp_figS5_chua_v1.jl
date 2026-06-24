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
    #size=(663,1000),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
    size=(400,620),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
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

optimal_coupling_regions_data_x_pos = DelimitedFiles.readdlm(DrWatson.datadir("chua","couplingregions_equidist__dist_0.20__cfrac_0.01__optimal_coupling__alpha_1.00_0.00_0.00.tsv"), '\t', Float64)
optimal_coupling_regions_data_y_pos = DelimitedFiles.readdlm(DrWatson.datadir("chua","couplingregions_equidist__dist_0.20__cfrac_0.01__optimal_coupling__alpha_0.00_1.00_0.00.tsv"), '\t', Float64)
optimal_coupling_regions_data_z_pos = DelimitedFiles.readdlm(DrWatson.datadir("chua","couplingregions_equidist__dist_0.20__cfrac_0.01__optimal_coupling__alpha_0.00_0.00_1.00.tsv"), '\t', Float64)

optimal_coupling_regions_data_x_neg = DelimitedFiles.readdlm(DrWatson.datadir("chua","couplingregions_equidist__dist_0.20__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00.tsv"), '\t', Float64)
optimal_coupling_regions_data_y_neg = DelimitedFiles.readdlm(DrWatson.datadir("chua","couplingregions_equidist__dist_0.20__cfrac_0.01__optimal_coupling__alpha_0.00_-1.00_0.00.tsv"), '\t', Float64)
optimal_coupling_regions_data_z_neg = DelimitedFiles.readdlm(DrWatson.datadir("chua","couplingregions_equidist__dist_0.20__cfrac_0.01__optimal_coupling__alpha_0.00_0.00_-1.00.tsv"), '\t', Float64)

default_color = TUD_darkblue
uncoupled_color = TUD_grey60
optimal_coupling_region_color = TUD_green1
optimal_coupling_region_timeseries_color = RGB(70,178,125)
coupling_region_example_1_color = TUD_orange1
coupling_region_example_2_color = TUD_orange2
coupling_region_example_3_color = TUD_red1

attractor_plot_range = 1:25000
number_of_optimal_coupling_regions = 15




#structure of figure

fig = Figure(
            #image size in Makie units converted with 300dpi from 183mm x 90mm (using Makies 0.75 pt_per_unit)
            #size = (ceil(183*25.4*300/0.75), ceil(90*25.4*300/0.75))
        )

#overall figure layout
group_fig = fig[1,1] = GridLayout()

panel_title = group_fig[1,1] = GridLayout()
panel_content = group_fig[2,1] = GridLayout()

panel_empty = panel_title[1,1] = GridLayout()
panel_pos = panel_title[1,2] = GridLayout()
panel_neg = panel_title[1,3] = GridLayout()

panel_x = panel_content[1,1] = GridLayout()
panel_y = panel_content[2,1] = GridLayout()
panel_z = panel_content[3,1] = GridLayout()

panel_a = panel_content[1,2] = GridLayout()
panel_b = panel_content[2,2] = GridLayout()
panel_c = panel_content[3,2] = GridLayout()
panel_d = panel_content[1,3] = GridLayout()
panel_e = panel_content[2,3] = GridLayout()
panel_f = panel_content[3,3] = GridLayout()

colsize!(panel_title,1,Auto(0.04))
colsize!(panel_title,2,Auto(0.48))
colsize!(panel_title,3,Auto(0.48))

colsize!(panel_content,1,Auto(0.04))
colsize!(panel_content,2,Auto(0.48))
colsize!(panel_content,3,Auto(0.48))

rowsize!(group_fig,1,Auto(0.04))
rowsize!(group_fig,2,Auto(0.96))

rowsize!(panel_content,1,Auto(0.32))
rowsize!(panel_content,2,Auto(0.32))
rowsize!(panel_content,3,Auto(0.32))



#individual panels

#panel title




Label(
    panel_pos[1, 1], 
    L"\text{coupling strength}\;\;\alpha = 1.0", 
    fontsize = 12,
    tellwidth = false, 
    tellheight = false
)

Label(
    panel_neg[1, 1], 
    L"\text{coupling strength}\;\;\alpha = -1.0", 
    fontsize = 12,
    tellwidth = false, 
    tellheight = false
)

Label(
    panel_x[1, 1], 
    #"x - coupling",
    #fontsize = 12,
    L"x\text{-coupling}",
    fontsize = 12,
    tellwidth = false, 
    tellheight = false,
    rotation = pi/2
)

Label(
    panel_y[1, 1], 
    L"y\text{-coupling}",
    fontsize = 12,
    tellwidth = false, 
    tellheight = false,
    rotation = pi/2
)

Label(
    panel_z[1, 1], 
    L"z\text{-coupling}",
    fontsize = 12,
    tellwidth = false, 
    tellheight = false,
    rotation = pi/2
)


#panel a: attractor plot optimal coupling regions

axis_panel_a = 
Axis3(
    panel_a[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.85 * pi,
    elevation = 0.15 * pi,
    limits = (-2.75,2.75,-2.75,2.75,-3.75,3.75),
    xticks = [-2,-1,0,1,2],
    xtickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    yticks = [-2,-1,0,1,2],
    ytickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    zticks = [-3,-2,-1,0,1,2,3],
    ztickformat = values -> [L"-3",L"-2",L"-1",L"0",L"1",L"2",L"3"],
    protrusions = (0,27,15,5),
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

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("chua","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_a, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_a, 
    optimal_coupling_regions_data_x_pos[1:number_of_optimal_coupling_regions,1],
    optimal_coupling_regions_data_x_pos[1:number_of_optimal_coupling_regions,2],
    optimal_coupling_regions_data_x_pos[1:number_of_optimal_coupling_regions,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data_x_pos[k,1], optimal_coupling_regions_data_x_pos[k,2], optimal_coupling_regions_data_x_pos[k,3]) for k in 1:number_of_optimal_coupling_regions]
meshscatter!(axis_panel_a,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data_x_pos[1:number_of_optimal_coupling_regions,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)

text!(axis_panel_a, 
    (0.02, 1.0), 
    text = L"\lambda \approx 0.514", #0.5144449439807204
    space = :relative, 
    align = (:left, :top), 
    color = :black,
    fontsize = 12
    )  

Label(
    panel_a[1, 1, TopLeft()], 
    "a",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -3*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel b: attractor plot optimal coupling regions

axis_panel_b = 
Axis3(
    panel_b[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.85 * pi,
    elevation = 0.15 * pi,
    limits = (-2.75,2.75,-2.75,2.75,-3.75,3.75),
    xticks = [-2,-1,0,1,2],
    xtickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    yticks = [-2,-1,0,1,2],
    ytickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    zticks = [-3,-2,-1,0,1,2,3],
    ztickformat = values -> [L"-3",L"-2",L"-1",L"0",L"1",L"2",L"3"],
    protrusions = (0,27,15,5),
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

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("chua","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_b, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_b, 
    optimal_coupling_regions_data_y_pos[1:number_of_optimal_coupling_regions,1],
    optimal_coupling_regions_data_y_pos[1:number_of_optimal_coupling_regions,2],
    optimal_coupling_regions_data_y_pos[1:number_of_optimal_coupling_regions,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data_y_pos[k,1], optimal_coupling_regions_data_y_pos[k,2], optimal_coupling_regions_data_y_pos[k,3]) for k in 1:number_of_optimal_coupling_regions]
meshscatter!(axis_panel_b,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data_y_pos[1:number_of_optimal_coupling_regions,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)

text!(axis_panel_b, 
    (0.02, 1.0), 
    text = L"\lambda \approx 0.455", #0.45505559365124926
    space = :relative, 
    align = (:left, :top), 
    color = :black,
    fontsize = 12
    )  

Label(
    panel_b[1, 1, TopLeft()], 
    "b",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -3*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel b: attractor plot optimal coupling regions

axis_panel_c = 
Axis3(
    panel_c[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.85 * pi,
    elevation = 0.15 * pi,
    limits = (-2.75,2.75,-2.75,2.75,-3.75,3.75),
    xticks = [-2,-1,0,1,2],
    xtickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    yticks = [-2,-1,0,1,2],
    ytickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    zticks = [-3,-2,-1,0,1,2,3],
    ztickformat = values -> [L"-3",L"-2",L"-1",L"0",L"1",L"2",L"3"],
    protrusions = (0,27,15,5),
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

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("chua","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_c, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_c, 
    optimal_coupling_regions_data_z_pos[1:number_of_optimal_coupling_regions,1],
    optimal_coupling_regions_data_z_pos[1:number_of_optimal_coupling_regions,2],
    optimal_coupling_regions_data_z_pos[1:number_of_optimal_coupling_regions,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data_z_pos[k,1], optimal_coupling_regions_data_z_pos[k,2], optimal_coupling_regions_data_z_pos[k,3]) for k in 1:number_of_optimal_coupling_regions]
meshscatter!(axis_panel_c,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data_z_pos[1:number_of_optimal_coupling_regions,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)

text!(axis_panel_c, 
    (0.02, 1.0), 
    text = L"\lambda \approx 0.437", #0.43741667061473477
    space = :relative, 
    align = (:left, :top), 
    color = :black,
    fontsize = 12
    )  

Label(
    panel_c[1, 1, TopLeft()], 
    "c",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -3*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)






#panel d: attractor plot optimal coupling regions

axis_panel_d = 
Axis3(
    panel_d[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.85 * pi,
    elevation = 0.15 * pi,
    limits = (-2.75,2.75,-2.75,2.75,-3.75,3.75),
    xticks = [-2,-1,0,1,2],
    xtickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    yticks = [-2,-1,0,1,2],
    ytickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    zticks = [-3,-2,-1,0,1,2,3],
    ztickformat = values -> [L"-3",L"-2",L"-1",L"0",L"1",L"2",L"3"],
    protrusions = (0,27,15,5),
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

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("chua","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_d, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_d, 
    optimal_coupling_regions_data_x_neg[1:number_of_optimal_coupling_regions,1],
    optimal_coupling_regions_data_x_neg[1:number_of_optimal_coupling_regions,2],
    optimal_coupling_regions_data_x_neg[1:number_of_optimal_coupling_regions,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data_x_neg[k,1], optimal_coupling_regions_data_x_neg[k,2], optimal_coupling_regions_data_x_neg[k,3]) for k in 1:number_of_optimal_coupling_regions]
meshscatter!(axis_panel_d,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data_x_neg[1:number_of_optimal_coupling_regions,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)

text!(axis_panel_d, 
    (0.02, 1.0), 
    text = L"\lambda \approx 0.531", #0.530670524256471
    space = :relative, 
    align = (:left, :top), 
    color = :black,
    fontsize = 12
    )  

Label(
    panel_d[1, 1, TopLeft()], 
    "d",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -3*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel e: attractor plot optimal coupling regions

axis_panel_e = 
Axis3(
    panel_e[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.85 * pi,
    elevation = 0.15 * pi,
    limits = (-2.75,2.75,-2.75,2.75,-3.75,3.75),
    xticks = [-2,-1,0,1,2],
    xtickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    yticks = [-2,-1,0,1,2],
    ytickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    zticks = [-3,-2,-1,0,1,2,3],
    ztickformat = values -> [L"-3",L"-2",L"-1",L"0",L"1",L"2",L"3"],
    protrusions = (0,27,15,5),
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

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("chua","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_e, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_e, 
    optimal_coupling_regions_data_y_neg[1:number_of_optimal_coupling_regions,1],
    optimal_coupling_regions_data_y_neg[1:number_of_optimal_coupling_regions,2],
    optimal_coupling_regions_data_y_neg[1:number_of_optimal_coupling_regions,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data_y_neg[k,1], optimal_coupling_regions_data_y_neg[k,2], optimal_coupling_regions_data_y_neg[k,3]) for k in 1:number_of_optimal_coupling_regions]
meshscatter!(axis_panel_e,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data_y_neg[1:number_of_optimal_coupling_regions,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)

text!(axis_panel_e, 
    (0.02, 1.0), 
    text = L"\lambda \approx 0.540", #0.5397005601082433
    space = :relative, 
    align = (:left, :top), 
    color = :black,
    fontsize = 12
    )  

Label(
    panel_e[1, 1, TopLeft()], 
    "e",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -3*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel f: attractor plot optimal coupling regions

axis_panel_f = 
Axis3(
    panel_f[1, 1],
    xlabel = L"x", 
    ylabel = L"y", 
    zlabel = L"z",
    aspect = :equal,
    azimuth = 1.85 * pi,
    elevation = 0.15 * pi,
    limits = (-2.75,2.75,-2.75,2.75,-3.75,3.75),
    xticks = [-2,-1,0,1,2],
    xtickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    yticks = [-2,-1,0,1,2],
    ytickformat = values -> [L"-2",L"-1",L"0",L"1",L"2"],
    zticks = [-3,-2,-1,0,1,2,3],
    ztickformat = values -> [L"-3",L"-2",L"-1",L"0",L"1",L"2",L"3"],
    protrusions = (0,27,15,5),
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

attractor_plot_data = DelimitedFiles.readdlm(DrWatson.datadir("chua","attractor_plot.tsv"), '\t', Float64)


lines!(axis_panel_f, attractor_plot_data[attractor_plot_range,2:4],
    linewidth = 0.33, 
    color = TUD_darkblue,
    alpha = 0.33,
    transparency = true
)

#plot optimal coupling sphere example
scatter!(
    axis_panel_f, 
    optimal_coupling_regions_data_z_neg[1:number_of_optimal_coupling_regions,1],
    optimal_coupling_regions_data_z_neg[1:number_of_optimal_coupling_regions,2],
    optimal_coupling_regions_data_z_neg[1:number_of_optimal_coupling_regions,3],
    color = optimal_coupling_region_color, 
    markersize = 5,      # Set a clear marker size
    strokewidth = 0.0,     # Remove stroke for cleaner appearance
    alpha = 1.0,
    transparency=true
)
spherecenter_positive = [Point(optimal_coupling_regions_data_z_neg[k,1], optimal_coupling_regions_data_z_neg[k,2], optimal_coupling_regions_data_z_neg[k,3]) for k in 1:number_of_optimal_coupling_regions]
meshscatter!(axis_panel_f,
    spherecenter_positive,
    marker = Sphere(Point(0.0,0.0,0.0),1),
    markersize = optimal_coupling_regions_data_z_neg[1:number_of_optimal_coupling_regions,4], 
    color = optimal_coupling_region_color,
    alpha = 0.25,
    transparency = true,
    rasterize = 10
)

text!(axis_panel_f, 
    (0.02, 1.0), 
    text = L"\lambda \approx 0.540", #0.5404266433534557
    space = :relative, 
    align = (:left, :top), 
    color = :black,
    fontsize = 12
    )  

Label(
    panel_f[1, 1, TopLeft()], 
    "f",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 5*unit_per_mm, -3*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




save(plotsdir("FIGS5_chua_collection_v1.pdf"), fig)
save(plotsdir("FIGS5_chua_collection_v1.png"), fig)


fig
