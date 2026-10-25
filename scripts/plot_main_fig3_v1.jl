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
    size=(400,300),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
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

# uncomment to rerun linearized comparison trajectory
#simulate_linearized_trajectory(roesslerODE_coupled_linear,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4], init_state[5], init_state[6]], transient_time = 0.0, sampling = true, output_filename = "linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 0.01, simulation_time = 10000.0)


# uncomment to rerun nonliear comparison trajectory
#simulate_nonlinear_trajectory(roesslerODE_coupled,"roessler",systemparameters_uncoupled; initial_state_index = -1, fixed_x0 = StaticArrays.@SVector[init_state[1], init_state[2], init_state[3]], delta_x0 = StaticArrays.@SVector[init_state[4]-init_state[1], init_state[5]-init_state[2], init_state[6]-init_state[3]], simulation_time = 1.0e6, transient_time = 0.0, sampling = true, output_filename = "nonlinear__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv", sampling_timestep = 1.0)


#structure of figure

fig = Figure(
            #image size in Makie units converted with 300dpi from 183mm x 90mm (using Makies 0.75 pt_per_unit)
            #size = (ceil(183*25.4*300/0.75), ceil(90*25.4*300/0.75))
        )

#overall figure layout
group_fig = fig[1,1] = GridLayout()

panel_a = group_fig[1,1] = GridLayout()
panel_b = group_fig[2,1] = GridLayout()
panel_c = group_fig[1,2] = GridLayout()
panel_d = group_fig[2,2] = GridLayout()




colsize!(group_fig,1,Auto(0.5))
colsize!(group_fig,2,Auto(0.5))

rowsize!(group_fig,1,Auto(0.5))
rowsize!(group_fig,2,Auto(0.5))

colgap!(group_fig, 4*unit_per_mm)



#individual panels



#panel a: coupled uncoupled comparison plot

synch_error_data_opt_full = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv"), '\t', Float64)
compute_stable_direction_uncoupled_linearized(
    roessler_backward_linear,
    "roessler", 
    "linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot.tsv",
    transient_time = 900.0
)
synch_error_data_opt_full_stabledirecton = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00_ncr_15__trajectory_plot__stabledirection.tsv"), '\t', Float64)


function is_coupled(x,y,z)    
    for k in 1:number_of_optimal_coupling_regions
        dist_sq = (x - optimal_coupling_regions_data[k,1])^2 + 
                (y - optimal_coupling_regions_data[k,2])^2 + 
                (z - optimal_coupling_regions_data[k,3])^2
        if dist_sq < optimal_coupling_regions_data[k,4]^2
            return true
        end
    end
    return false
end

synch_error_data_opt_full_active_coupling_length = size(synch_error_data_opt_full,1)
synch_error_data_opt_full_stabledirecton_length = size(synch_error_data_opt_full_stabledirecton,1)
#synch_error_data_opt_full_active_coupling_length = 10000
synch_error_data_opt_full_active_coupling = [ 
    #is_in_coupling_region("roessler", synch_error_data_opt_full[k,2:4], "couplingregions_equidist__dist_0.50__cfrac_0.01__optimal_coupling__alpha_-1.00_0.00_0.00.tsv", number_of_optimal_coupling_regions) 
    is_coupled(synch_error_data_opt_full[k,2], synch_error_data_opt_full[k,3], synch_error_data_opt_full[k,4])
    for k in 1:synch_error_data_opt_full_active_coupling_length
]

t_full = synch_error_data_opt_full[1:synch_error_data_opt_full_active_coupling_length,1]
deltax_full = exp.(synch_error_data_opt_full[1:synch_error_data_opt_full_active_coupling_length,8])
deltax_highlight = [synch_error_data_opt_full_active_coupling[k] ? deltax_full[k] : NaN for k in 1:synch_error_data_opt_full_active_coupling_length]


deltax_alignment_full = [
    abs(LinearAlgebra.dot(synch_error_data_opt_full[k, 5:7], synch_error_data_opt_full_stabledirecton[k, 2:4])) / LinearAlgebra.norm(synch_error_data_opt_full[k, 5:7]) / LinearAlgebra.norm(synch_error_data_opt_full_stabledirecton[k, 2:4]) 
    for k in 1:synch_error_data_opt_full_stabledirecton_length
]
deltax_alignment_highlight = [synch_error_data_opt_full_active_coupling[k] ? deltax_alignment_full[k] : NaN for k in 1:synch_error_data_opt_full_stabledirecton_length]

mask_full = .!isnan.(deltax_highlight[1:synch_error_data_opt_full_active_coupling_length])
# Pad the mask with false to detect edges at the boundaries
padded_full = [false; mask_full; false]
# A block starts if current is true and previous was false
coupling_indices = findall(i -> padded_full[i] && !padded_full[i-1], 2:length(padded_full))


mask = .!isnan.(deltax_alignment_highlight[1:synch_error_data_opt_full_stabledirecton_length])
# Pad the mask with false to detect edges at the boundaries
padded = [false; mask; false]
# A block starts if current is true and previous was false
starts = findall(i -> padded[i] && !padded[i-1], 2:length(padded))
# A block ends if current is true and next is false
ends = findall(i -> padded[i] && !padded[i+1], 2:length(padded))
datablocks = collect(zip(starts, ends))




synch_error_data_uncoupled = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv"), '\t', Float64)
compute_stable_direction_uncoupled_linearized(
    roessler_backward_linear,
    "roessler", 
    "linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison.tsv",
    transient_time = 900.0
)
synch_error_data_uncoupled_stabledirecton = DelimitedFiles.readdlm(DrWatson.datadir("roessler","linearized__xcoupled_trajectory_plot__alpha_-1.00_0.00_0.00__uncoupled_comparison__stabledirection.tsv"), '\t', Float64)

deltax_alignment_full_uncoupled = [
    abs(LinearAlgebra.dot(synch_error_data_uncoupled[k, 5:7], synch_error_data_uncoupled_stabledirecton[k, 2:4])) / LinearAlgebra.norm(synch_error_data_uncoupled[k, 5:7]) / LinearAlgebra.norm(synch_error_data_uncoupled_stabledirecton[k, 2:4]) 
    for k in 1:size(synch_error_data_uncoupled_stabledirecton,1)
]




axis_panel_a = 
Axis(
    panel_a[1, 1],
    ylabel = L"\text{synch. error}\;\left|\left| \Delta \mathbf{x} \, \right|\right|",
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
    axis_panel_a,
    synch_error_data_uncoupled[1:10000,1],exp.(synch_error_data_uncoupled[1:10000,8]),
    color = uncoupled_color,
    linewidth = 1.5
)

lines!(
    axis_panel_a,
    t_full[1:10000], deltax_full[1:10000],
    color = default_color,
    linewidth = 1.5
)


scatter!(
    axis_panel_a,
    t_full[1:10000], deltax_highlight[1:10000],
    color = optimal_coupling_region_timeseries_color,
    markersize = 3.0
)

for (s,e) in datablocks[1:100]
    vspan!(
        axis_panel_a, 
        t_full[s], 
        t_full[e], 
        color = optimal_coupling_region_timeseries_color,
        alpha = 0.2
    )
end


text!(axis_panel_a, 
    (7, 1e1), 
    text = L"\textbf{A}", 
    #space = :relative,
    align = (:center, :center),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_a,
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


text!(axis_panel_a, 
    (5, 4e-2), 
    text = L"\textbf{B,C}", 
    #space = :relative,
    align = (:center, :center),
    color = optimal_coupling_region_color,
    fontsize = 14
    ) 

arrows2d!(axis_panel_a,
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

arrows2d!(axis_panel_a,
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
    panel_a[1, 1, TopLeft()], 
    "a",
    fontsize = 12,
    font = :bold,
    padding = (-4*unit_per_mm, 8*unit_per_mm, -3*unit_per_mm, -3*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel b: alignment with stable direction

axis_panel_b = 
Axis(
    panel_b[1, 1],
    #ylabel = L"\left|\Delta x\,\right|",
    ylabel = L"\left|\Delta \mathbf{x} \cdot \Delta \mathbf{x}_\mathrm{stable}\,\right|",
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
    axis_panel_b,
    t_full[1:10000], deltax_alignment_full_uncoupled[1:10000],
    color = uncoupled_color,
    linewidth = 1.5
)

lines!(
    axis_panel_b,
    t_full[1:10000], deltax_alignment_full[1:10000],
    color = default_color,
    linewidth = 1.5
)

scatter!(
    axis_panel_b,
    t_full[1:10000], deltax_alignment_highlight[1:10000],
    color = optimal_coupling_region_timeseries_color,
    markersize = 3.0
)


for (s,e) in datablocks[1:100]
    vspan!(
        axis_panel_b, 
        t_full[s], 
        t_full[e], 
        color = optimal_coupling_region_timeseries_color,
        alpha = 0.2
    )
end


Label(
    panel_b[1, 1, TopLeft()], 
    "b",
    fontsize = 12,
    font = :bold,
    padding = (-4*unit_per_mm, 8*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom
)




#panel c: instantaneous contributions to synchronization error

roessler_jacobian(x,y,z) = [0 -1 -1; 1 0.2 0; z 0 (x-5.7)]

error_growth_intrinsic_contribution_full = [ 
        let x = state[2], y = state[3], z = state[4], dx = state[5], dy = state[6], dz = state[7]
            LinearAlgebra.dot(roessler_jacobian(x,y,z)*[dx, dy, dz],[dx, dy, dz]) 
        end
        for state in eachrow(synch_error_data_opt_full)
    ]
error_growth_coupling_contribution_full = [ 
        let x = state[2], y = state[3], z = state[4], dx = state[5], dy = state[6], dz = state[7]
            is_coupled(x,y,z) ? dx*dx : 0.0
        end
        for state in eachrow(synch_error_data_opt_full)
    ]

error_growth_intrinsic_contribution_uncoupled_comparison = [ 
        let x = state[2], y = state[3], z = state[4], dx = state[5], dy = state[6], dz = state[7]
            LinearAlgebra.dot(roessler_jacobian(x,y,z)*[dx, dy, dz],[dx, dy, dz]) 
        end
        for state in eachrow(synch_error_data_uncoupled)
    ]


axis_panel_c = 
Axis(
    panel_c[1, 1],
    #ylabel = L"\left|\Delta x\,\right|",
    ylabel = L"\text{expansion rate}\;\lambda^\mathrm{inst}(t)",
    yscale = Makie.Symlog10(0.1),
    limits = (0,22.5,-20,20),
    xticks = [0,5,10,15,20],
    xtickformat = values -> [L"0",L"5",L"10",L"15",L"20"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 2.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [-10, -1, -0.1, 0, 0.1, 1, 10],
    ytickformat = values -> [L"-10", L"-1", L"-0.1", L"0", L"0.1", L"1", L"10"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

lines!(
    axis_panel_c,
    t_full[1:10000], error_growth_intrinsic_contribution_uncoupled_comparison[1:10000],
    color = uncoupled_color,
    linewidth = 1.5
)

lines!(
    axis_panel_c,
    t_full[1:10000], error_growth_intrinsic_contribution_full[1:10000],
    color = default_color,
    linewidth = 1.5
)

lines!(
    axis_panel_c,
    t_full[1:10000], error_growth_coupling_contribution_full[1:10000],
    color = optimal_coupling_region_color,
    linewidth = 1.5
)


for (s,e) in datablocks[1:100]
    vspan!(
        axis_panel_c, 
        t_full[s], 
        t_full[e], 
        color = optimal_coupling_region_timeseries_color,
        alpha = 0.2
    )
end


Label(
    panel_c[1, 1, TopLeft()], 
    "c",
    fontsize = 12,
    font = :bold,
    padding = (-8*unit_per_mm, 10*unit_per_mm, -3*unit_per_mm, -3*unit_per_mm),
    halign = :right,
    valign = :bottom
)



error_growth_total_full = error_growth_intrinsic_contribution_full + error_growth_coupling_contribution_full

tau = 2.0
time_series_dt = 0.01
tau_index_shift = round(Int, tau/time_series_dt)
finite_time_error_growth_total_full = [ (synch_error_data_opt_full[i + tau_index_shift, 8] - synch_error_data_opt_full[i, 8]) / tau for i in 1:(size(synch_error_data_opt_full,1) - tau_index_shift) ]
finite_time_error_growth_total_uncoupled_comparison= [ (synch_error_data_uncoupled[i + tau_index_shift, 8] - synch_error_data_uncoupled[i, 8]) / tau for i in 1:(size(synch_error_data_uncoupled,1) - tau_index_shift) ]

#avoid coupling events at the end of the trajectory where no finite time expansion rates can be computed
plot_coupling_indices = coupling_indices[1:end-10]

function is_in_area_A_region(x,y,z)    
    for k in 1:number_of_optimal_coupling_regions
        if(optimal_coupling_regions_data[k,1] < 4.0)
            continue
        end
        dist_sq = (x - optimal_coupling_regions_data[k,1])^2 + 
                (y - optimal_coupling_regions_data[k,2])^2 + 
                (z - optimal_coupling_regions_data[k,3])^2
        if dist_sq < optimal_coupling_regions_data[k,4]^2
            return true
        end
    end
    return false
end
area_A_coupling_indices = [i for i in plot_coupling_indices if is_in_area_A_region(synch_error_data_opt_full[i,2], synch_error_data_opt_full[i,3], synch_error_data_opt_full[i,4])]

axis_panel_d = 
Axis(
    panel_d[1, 1],
    #ylabel = L"\left|\Delta x\,\right|",
    ylabel = L"\text{expansion rate}\;\lambda^\mathrm{FT}_{\tau=2}(t)",
    xscale = Makie.Symlog10(0.1),
    yscale = Makie.Symlog10(0.1),
    limits = (-20,20,-20,20),
    xticks = [-10, -1, -0.1, 0, 0.1, 1, 10],
    xtickformat = values -> [L"-10", L"-1", L"-0.1", L"0", L"0.1", L"1", L"10"],
    xlabel = L"\text{expansion rate}\;\lambda^\mathrm{inst}(t)",
    xlabelpadding = -1,
    ylabelpadding = 2.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [-10, -1, -0.1, 0, 0.1, 1, 10],
    ytickformat = values -> [L"-10", L"-1", L"-0.1", L"0", L"0.1", L"1", L"10"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

scatter!(
    axis_panel_d,
    error_growth_intrinsic_contribution_uncoupled_comparison[1:50:size(finite_time_error_growth_total_uncoupled_comparison,1)], finite_time_error_growth_total_uncoupled_comparison[1:50:end],
    color = uncoupled_color,
    markersize = 1.5
)

scatter!(
    axis_panel_d,
    error_growth_intrinsic_contribution_full[1:50:size(finite_time_error_growth_total_full,1)], finite_time_error_growth_total_full[1:50:end],
    color = default_color,
    markersize = 1
)


scatter!(
    axis_panel_d,
    [error_growth_total_full[i+1] for i in plot_coupling_indices], [finite_time_error_growth_total_full[i+1] for i in plot_coupling_indices],
    color = optimal_coupling_region_color,
    markersize = 2
)

scatter!(
    axis_panel_d,
    [error_growth_total_full[i+1] for i in area_A_coupling_indices], [finite_time_error_growth_total_full[i+1] for i in area_A_coupling_indices],
    color = optimal_coupling_region_timeseries_color,
    markersize = 2
)


Label(
    panel_d[1, 1, TopLeft()], 
    "d",
    fontsize = 12,
    font = :bold,
    padding = (-8*unit_per_mm, 10*unit_per_mm, -2*unit_per_mm, -2*unit_per_mm),
    halign = :right,
    valign = :bottom
)




save(plotsdir("FIG3_optimal_coupling_alpha_-1.00_mechanism_v1.pdf"), fig)
save(plotsdir("FIG3_optimal_coupling_alpha_-1.00_mechanism_v1.png"), fig)


fig
