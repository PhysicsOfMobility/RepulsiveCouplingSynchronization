using DrWatson
@quickactivate :ChaosSynchronization

import DelimitedFiles
import StaticArrays 
import DifferentialEquations 
import DynamicalSystems

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
    #size=(663,300),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
    size=(331,300),     #converted to 90 pt per inch (seems to more correctly reproduce size in LaTeX)
    figure_padding=1*unit_per_mm,  #1mm padding
    fontsize = 10,  #in pt
    fonts = (;
        regular = "Noto Sans", 
        bold = "Noto Sans Bold", 
        italic = "Noto Sans Italic", 
        bold_italic = "Noto Sans Bold Italic"
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


#TUD_theme = merge(TUD_theme, theme_latexfonts())

set_theme!(TUD_theme)




default_color = TUD_darkblue
uncoupled_color = TUD_grey60
optimal_coupling_region_color = TUD_green1
optimal_coupling_region_timeseries_color = RGB(70,178,125)
optimal_coupling_region_timeseries_color_lighter = RGB(142, 211, 177)
coupling_region_example_1_color = TUD_orange1
coupling_region_example_2_color = TUD_orange2
coupling_region_example_3_color = TUD_red1



#set up all data required for later figures

struct LimitCycleParameters{Float64}
    couplingstrength::Float64
    couplingcenter::Float64
    couplingregionssize::Float64
end

function limitcycleODE_coupled(x::StaticArrays.SVector{4, Float64}, p::LimitCycleParameters, t::Float64) 
    r1 = sqrt(x[1]^2 + x[2]^2)
    r2 = sqrt(x[3]^2 + x[4]^2)

    dx1 = x[1]/r1 - x[1] - x[2]*r1
    dx2 = x[2]/r1 - x[2] + x[1]*r1

    dx3 = x[3]/r2 - x[3] - x[4]*r2 
    dx4 = x[4]/r2 - x[4] + x[3]*r2 

    if abs(atan(x[2], x[1]) - p.couplingcenter) <= p.couplingregionssize || abs(atan(x[2], x[1]) - p.couplingcenter + 2*pi) <= p.couplingregionssize || abs(atan(x[2], x[1]) - p.couplingcenter - 2*pi) <= p.couplingregionssize
        dx3 -= p.couplingstrength * (x[3] - x[1])
    end

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4]
end

function limitcycleODE_coupled_linear(x::StaticArrays.SVector{5, Float64}, p::LimitCycleParameters, t::Float64) 
    r = sqrt(x[1]^2 + x[2]^2)

    dx1 = x[1]/r - x[1] - x[2]*r 
    dx2 = x[2]/r - x[2] + x[1]*r 

    dx3 = (-x[1]^2-x[1]*x[2]) * x[3] + (-1-x[1]*x[2]-x[2]^2) * x[4]
    dx4 = (1+x[1]^2-x[1]*x[2]) * x[3] + (x[1]*x[2]-x[2]^2) * x[4]

    if abs(atan(x[2], x[1]) - p.couplingcenter) <= p.couplingregionssize || abs(atan(x[2], x[1]) - p.couplingcenter + 2*pi) <= p.couplingregionssize || abs(atan(x[2], x[1]) - p.couplingcenter - 2*pi) <= p.couplingregionssize
        dx3 -= p.couplingstrength * x[3]
    end

    dx5 = dx3 * x[3] + dx4 * x[4]

    # subtract radial component
    dx3 = dx3 - x[3] * dx5
    dx4 = dx4 - x[4] * dx5

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5]
end

# Out-of-place (OOP) callback for immutable StaticArrays
function normalizecallback_affect!(integrator)
    u = integrator.u
    delta = sqrt(u[3]^2 + u[4]^2)
    v3 = u[3] / delta
    v4 = u[4] / delta
    # Reconstruct the SVector
    integrator.u = StaticArrays.@SVector[u[1], u[2], v3, v4, u[5]]
end
init_normalize = (c,u,t,integrator) -> normalizecallback_affect!(integrator) # Initialize state to be normalized
linearizedODE_normalizecallback = DifferentialEquations.DiscreteCallback((u,t,integrator) -> true, normalizecallback_affect!, initialize=init_normalize, save_positions=(false,false))


absolute_tolerance = 1e-6
relative_tolerance = 1e-4
max_integration_timestep = 1e-2

transient_time = 100*2*pi
simulation_time = 1000*2*pi

coupling_angle_range = 0:1:360.5

coupling_strength = -1.0



#=
layp_pos_coupling = zeros(Float64, (size(coupling_angle_range,1),2))

for c in 1:size(coupling_angle_range,1)
    p = LimitCycleParameters(
        -coupling_strength, 
        2*pi*coupling_angle_range[c]/360.0,
        pi/32.0
    )

    x0 = StaticArrays.SVector{5,Float64}(1.0, 0.0, sqrt(0.5), sqrt(0.5), 0.0)

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            limitcycleODE_coupled_linear, 
            x0,
            p,
            diffeq = (
                callback = linearizedODE_normalizecallback,
                save_everystep = false, 
                save_start = false,
                dtmax = max_integration_timestep,
                abstol = absolute_tolerance,
                reltol = relative_tolerance,
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)
    x_transient = DynamicalSystems.current_state(ds)
    x0_reinit = StaticArrays.setindex(x_transient, 0.0, 5) # Use setindex
    DynamicalSystems.reinit!(ds, x0_reinit, t0 = 0.0)

    DynamicalSystems.step!(ds, simulation_time, true)

    layp_pos_coupling[c,:] = [2*pi*coupling_angle_range[c]/360.0, DynamicalSystems.current_state(ds)[5]/DynamicalSystems.current_time(ds)]
end
=#


layp_neg_coupling = zeros(Float64, (size(coupling_angle_range,1),2))

for c in 1:size(coupling_angle_range,1)
    p = LimitCycleParameters(
        coupling_strength, 
        2*pi*coupling_angle_range[c]/360.0,
        pi/32.0
    )

    x0 = StaticArrays.SVector{5,Float64}(1.0, 0.0, sqrt(0.5), sqrt(0.5), 0.0)

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            limitcycleODE_coupled_linear, 
            x0,
            p,
            diffeq = (
                callback = linearizedODE_normalizecallback,
                save_everystep = false, 
                save_start = false,
                dtmax = max_integration_timestep,
                abstol = absolute_tolerance,
                reltol = relative_tolerance,
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)
    x_transient = DynamicalSystems.current_state(ds)
    x0_reinit = StaticArrays.setindex(x_transient, 0.0, 5) # Use setindex
    DynamicalSystems.reinit!(ds, x0_reinit, t0 = 0.0)

    DynamicalSystems.step!(ds, simulation_time, true)

    layp_neg_coupling[c,:] = [2*pi*coupling_angle_range[c]/360.0, DynamicalSystems.current_state(ds)[5]/DynamicalSystems.current_time(ds)]
end





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



#colsize!(group_fig, 1, Aspect(2.0, 1.0))
colsize!(group_fig,1,Auto(0.5))
colsize!(group_fig,2,Auto(0.5))

rowsize!(group_fig,1,Auto(0.5))
rowsize!(group_fig,2,Auto(0.5))

colgap!(group_fig, 2*unit_per_mm)
rowgap!(group_fig, 2*unit_per_mm)

#individual panels



#panel a: limit cycle illustration

axis_panel_a =
Axis(
    panel_a[1,1],
    limits = (-1.5,2,-1.5,2),
    aspect = 1,
    width = 175,   # Fixed size in pixels
    height = 175,  # Keeps it a perfect circle
    tellwidth = false, 
    tellheight = false,
    halign = :right
)
hidedecorations!(axis_panel_a)
hidespines!(axis_panel_a)

poly!(
    axis_panel_a, 
    Circle(Point2f(0, 0), 1), 
    color = :white, 
    strokecolor = :black, 
    #linestyle = (:dot,:dense),
    strokewidth = 2
)

poly!(
    axis_panel_a, 
    [Point2f(0.0,0.0), [Point2f(1.5 * cos(a), 1.5 * sin(a)) for a in 0:0.01:pi/4]... , Point2f(0.0,0.0)], 
    color = optimal_coupling_region_timeseries_color,
    alpha = 0.2,
    strokewidth = 0
)

poly!(
    axis_panel_a, 
    [Point2f(0.0,0.0), [Point2f(1.2 * cos(a), 1.2 * sin(a)) for a in pi:0.01:pi+pi/4]... , Point2f(0.0,0.0)], 
    color = optimal_coupling_region_timeseries_color_lighter,
    alpha = 0.2,
    strokewidth = 0
)

arrows2d!(
    axis_panel_a, 
    [0.0], 
    [0.0], 
    [1.75], 
    [0.0], 
    color = TUD_black,
    shaftwidth = 2,
    tipwidth = 6
)

text!(
    axis_panel_a, 
    1.55, 
    -0.25; 
    text = L"r", 
    align = (:center, :center), 
    fontsize = 12
)

lines!(
    axis_panel_a,
    [0,sqrt(0.75)],[0,sqrt(0.25)],
    color = default_color,
    linewidth = 2
)

text!(
    axis_panel_a, 
    0.7, 
    0.25; 
    text = L"\theta", 
    align = (:center, :center), 
    fontsize = 12
)



scatter!(
    axis_panel_a,
    [sqrt(0.75)],[sqrt(0.25)],
    color = default_color,
    markersize = 10
)

poly!(
    axis_panel_a, 
    Circle(Point2f(-sqrt(0.65), sqrt(0.35)), 0.075), 
    color = :white, 
    strokecolor = default_color,
    strokewidth = 1.5
)

lines!(
    axis_panel_a,
    [sqrt(0.75),sqrt(0.75)],[0.0,1.0],
    color = TUD_grey40,
    linewidth = 1.5,
    linestyle = (:dash, :dense)
)



scatter!(
    axis_panel_a,
    [sqrt(0.95)],[sqrt(0.05)],
    color = TUD_orange1,
    markersize = 10
)

arrows2d!(
    axis_panel_a, 
    [sqrt(0.75)], 
    [sqrt(0.05)], 
    [0.5], 
    [0], 
    color = TUD_orange1,
    shaftwidth = 2.5, 
    tipwidth = 9
)

scatter!(
    axis_panel_a,
    [sqrt(0.45)],[sqrt(0.55)],
    color = optimal_coupling_region_timeseries_color,
    markersize = 10
)

arrows2d!(
    axis_panel_a, 
    [sqrt(0.75)], 
    [sqrt(0.55)], 
    [-0.6], 
    [0], 
    color = optimal_coupling_region_timeseries_color,
    shaftwidth = 2.5, 
    tipwidth = 9
)

slow_r_list = range(0.75,step = (0.85-0.75)/100.0, length = 101)
slow_phi_list = range(1.3,step = (3*pi/4-1.3)/100.0, length = 101)

slow_x_list = [r*cos(p) for (r,p) in zip(slow_r_list,slow_phi_list)]
slow_y_list = [r*sin(p) for (r,p) in zip(slow_r_list,slow_phi_list)]

lines!(
    axis_panel_a,
    slow_x_list,slow_y_list,
    color = optimal_coupling_region_timeseries_color,
    linewidth = 2
)

arrows2d!(
    axis_panel_a, 
    [slow_x_list[end]+0.025], 
    [slow_y_list[end]+0.025], 
    [-0.1], 
    [-0.1], 
    tiplength = 30,
    tipwidth = 40,
    color = optimal_coupling_region_timeseries_color,
)

text!(
    axis_panel_a, 
    -0.1, 
    0.5; 
    text = "slow", 
    align = (:center, :center), 
    color = optimal_coupling_region_timeseries_color,
    fontsize = 12
)

fast_r_list = range(1.5,step = (1.15-1.5)/100.0, length = 101)
fast_phi_list = range(0.15,step = (3*pi/4-0.15)/100.0, length = 101)

fast_x_list = [r*cos(p) for (r,p) in zip(fast_r_list,fast_phi_list)]
fast_y_list = [r*sin(p) for (r,p) in zip(fast_r_list,fast_phi_list)]

lines!(
    axis_panel_a,
    fast_x_list,fast_y_list,
    color = TUD_orange1,
    linewidth = 2
)

arrows2d!(
    axis_panel_a, 
    [fast_x_list[end]+0.025], 
    [fast_y_list[end]+0.025], 
    [-0.1], 
    [-0.1], 
    tiplength = 30,
    tipwidth = 40,
    color = TUD_orange1,
)

text!(
    axis_panel_a, 
    1.2, 
    1.1; 
    text = "fast", 
    align = (:center, :center), 
    color = TUD_orange1,
    rotation = -0.225*pi,
    fontsize = 12
)


Label(
    panel_a[1, 1, TopLeft()], 
    "a",
    fontsize = 12,
    #color = TUD_grey40,
    font = :bold,
    padding = (-8*unit_per_mm, 11*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)






axis_panel_b = 
Axis(
    panel_b[1, 1],
    ylabel = L"\text{max. Lyap. exp.}\;\;\lambda",
    limits = (0,2*pi, -0.015,0.045),
    xticks = [0,pi/2,pi,3*pi/2,2*pi,],
    xtickformat = values -> [L"0",L"\pi/2",L"\pi",L"3\,\pi/2",L"2\,\pi"],
    xlabel = L"\text{coupling center}\;\theta_c",
    xlabelpadding = -1,
    ylabelpadding = -1.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    yticks = [-0.01,0,0.01,0.02,0.03,0.04],
    ytickformat = values -> [L"-0.01",L"0",L"0.01",L"0.02",L"0.03",L"0.04"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

lines!(
    axis_panel_b,
    [0,2*pi],[0.0,0.0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

lines!(
    axis_panel_b, 
    layp_neg_coupling,
    color = default_color,
    linewidth = 2
)

scatter!(
    axis_panel_b,
    [pi/8],[-0.007010937880372775],
    color = default_color,
    markersize = 10
)

vspan!(
    axis_panel_b, 
    0, 
    pi/4, 
    color = optimal_coupling_region_timeseries_color,
    alpha = 0.2
)

vspan!(
    axis_panel_b, 
    pi, 
    5*pi/4, 
    color = optimal_coupling_region_timeseries_color_lighter,
    alpha = 0.2
)

#text!(
#    axis_panel_b, 
#    4.8, 
#    0.0415; 
#    text = L"\alpha = -1.0", 
#    align = (:center, :center), 
#    fontsize = 12
#)

Label(
    panel_b[1, 1, TopLeft()], 
    "b",
    fontsize = 12,
    font = :bold,
    padding = (-8*unit_per_mm, 11*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)





transient_time = 100*pi
simulation_time = 20*pi

sampling_timestep = 0.01



p = LimitCycleParameters(
    coupling_strength, 
    pi/8,
    pi/32.0
)

x0 = StaticArrays.SVector{5,Float64}(1.0, 0.0, sqrt(0.5), sqrt(0.5), 0.0)

ds = 
    DynamicalSystems.ContinuousDynamicalSystem(
        limitcycleODE_coupled_linear, 
        x0,
        p,
        diffeq = (
            callback = linearizedODE_normalizecallback,
            save_everystep = false, 
            save_start = false,
            dtmax = max_integration_timestep,
            abstol = absolute_tolerance,
            reltol = relative_tolerance,
        )
    )

DynamicalSystems.step!(ds, transient_time, true)
x_transient = DynamicalSystems.current_state(ds)
x0_reinit = StaticArrays.setindex(x_transient, 0.0, 5) # Use setindex
DynamicalSystems.reinit!(ds, x0_reinit, t0 = 0.0)

linearized_example_u, linearized_example_t = DynamicalSystems.trajectory(ds, simulation_time, Δt = sampling_timestep)

linearized_example_active_coupling = [ 
    pi/8 - pi/32 <= atan(linearized_example_u[k,2],linearized_example_u[k,1]) <= pi/8 + pi/32
    for k in 1:size(linearized_example_u,1)
]

t_full = linearized_example_t
deltax_full = exp.(linearized_example_u[:,5])
deltax_highlight = [linearized_example_active_coupling[k] ? deltax_full[k] : NaN for k in 1:size(linearized_example_u,1)]


mask = .!isnan.(deltax_highlight)
# Pad the mask with false to detect edges at the boundaries
padded = [false; mask; false]
# A block starts if current is true and previous was false
starts = findall(i -> padded[i] && !padded[i-1], 2:length(padded))
# A block ends if current is true and next is false
ends = findall(i -> padded[i] && !padded[i+1], 2:length(padded))
datablocks = collect(zip(starts, ends))




axis_panel_c = 
Axis(
    panel_c[1, 1],
    ylabel = L"\text{synch. error}\;\left|\left| \Delta \mathbf{x} \, \right|\right|",
    yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0, 11*pi, 7.5e-1, 1.05e0),
    xticks = [0,2*pi, 4*pi, 6*pi, 8*pi, 10*pi],
    xtickformat = values -> [L"0",L"2\,\pi",L"4\,\pi",L"6\,\pi",L"8\,\pi",L"10\,\pi"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 4.0,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [0.8,0.9,1.0],
    ytickformat = values -> [L"0.8",L"0.9",L"1.0"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

lines!(
    axis_panel_c,
    [0.0,11*pi],[0.0,0.0],
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


lines!(
    axis_panel_c,
    [pi,9*pi],[0.935, 0.935*exp(-0.007007683544933368*8*pi)],
    color = optimal_coupling_region_color,
    linewidth = 1
)

text!(
    axis_panel_c,
    (2*pi,0.875),
    text = L"\lambda \approx -0.007",  #-0.007007683544933368
    rotation = -0.21*pi,
    color = optimal_coupling_region_color
)

text!(
    axis_panel_c, 
    8.4*pi, 
    1.025; 
    text = L"\alpha = -1.0", 
    align = (:center, :center), 
    fontsize = 12
)

Label(
    panel_c[1, 1, TopLeft()], 
    "c",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 8*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)







transient_time = 100*pi
simulation_time = 100*pi

sampling_timestep = 0.01



p = LimitCycleParameters(
    coupling_strength, 
    pi/8,
    pi/32.0
)

x0 = StaticArrays.SVector{4,Float64}(1.0, 0.0, cos(-0.1*pi), sin(-0.1*pi))
x0 = StaticArrays.SVector{4,Float64}(1.0, 0.0, 1e-6, 1e-6)

ds = 
    DynamicalSystems.ContinuousDynamicalSystem(
        limitcycleODE_coupled, 
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

nonlinear_example_u, nonlinear_example_t = DynamicalSystems.trajectory(ds, simulation_time, Δt = sampling_timestep)

nonlinear_example_active_coupling = [ 
    pi/8 - pi/32 <= atan(nonlinear_example_u[k,2],nonlinear_example_u[k,1]) <= pi/8 + pi/32
    for k in 1:size(nonlinear_example_u,1)
]

t_full = nonlinear_example_t
deltax_full = sqrt.( (nonlinear_example_u[:,3].-nonlinear_example_u[:,1]).^2 .+ (nonlinear_example_u[:,4].-nonlinear_example_u[:,2]).^2 )
deltax_highlight = [nonlinear_example_active_coupling[k] ? deltax_full[k] : NaN for k in 1:size(nonlinear_example_u,1)]


mask = .!isnan.(deltax_highlight)
# Pad the mask with false to detect edges at the boundaries
padded = [false; mask; false]
# A block starts if current is true and previous was false
starts = findall(i -> padded[i] && !padded[i-1], 2:length(padded))
# A block ends if current is true and next is false
ends = findall(i -> padded[i] && !padded[i+1], 2:length(padded))
datablocks = collect(zip(starts, ends))




axis_panel_d = 
Axis(
    panel_d[1, 1],
    ylabel = L"\text{synch. error}\;\left|\left| \Delta \mathbf{x} \, \right|\right|",
    #yscale = log10, # Logarithmic y-axis
    #limits = (0,1e3 ,1e-1, 1e6),
    limits = (0, 15*pi, 0, 1),
    xticks = [0, 4*pi, 8*pi, 12*pi],
    xtickformat = values -> [L"0",L"4\,\pi",L"8\,\pi",L"12\,\pi"],
    xlabel = L"\text{time}\;t",
    xlabelpadding = -1,
    ylabelpadding = 4.0,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    #yticks = [1e0,1e2,1e4,1e6],
    #ytickformat = values -> [L"10^0",L"10^2",L"10^4",L"10^6"],
    yticks = [0,0.2,0.4,0.6,0.8,1.0],
    ytickformat = values -> [L"0",L"0.2",L"0.4",L"0.6",L"0.8",L"1.0"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

lines!(
    axis_panel_d,
    [0.0,11*pi],[0.0,0.0],
    color = uncoupled_color,
    linewidth = 1.5
)

lines!(
    axis_panel_d,
    t_full, deltax_full,
    color = default_color,
    linewidth = 1.5
)

scatter!(
    axis_panel_d,
    t_full, deltax_highlight,
    color = optimal_coupling_region_timeseries_color,
    markersize = 3.0
)

Label(
    panel_d[1, 1, TopLeft()], 
    "d",
    fontsize = 12,
    font = :bold,
    padding = (-5*unit_per_mm, 8*unit_per_mm, -4*unit_per_mm, -6*unit_per_mm),
    halign = :right,
    valign = :bottom
)




transient_time = 100*2*pi
simulation_time = 1000*2*pi

coupling_strength_range = -20:0.1:20

lyap_scan = zeros(Float64, (size(coupling_strength_range,1),2))

for c in 1:size(coupling_strength_range,1)
    p = LimitCycleParameters(
        coupling_strength_range[c], 
        pi/8.0,
        pi/32.0
    )

    x0 = StaticArrays.SVector{5,Float64}(1.0, 0.0, sqrt(0.5), sqrt(0.5), 0.0)

    ds = 
        DynamicalSystems.ContinuousDynamicalSystem(
            limitcycleODE_coupled_linear, 
            x0,
            p,
            diffeq = (
                callback = linearizedODE_normalizecallback,
                save_everystep = false, 
                save_start = false,
                dtmax = max_integration_timestep,
                abstol = absolute_tolerance,
                reltol = relative_tolerance,
            )
        )

    DynamicalSystems.step!(ds, transient_time, true)
    x_transient = DynamicalSystems.current_state(ds)
    x0_reinit = StaticArrays.setindex(x_transient, 0.0, 5) # Use setindex
    DynamicalSystems.reinit!(ds, x0_reinit, t0 = 0.0)

    DynamicalSystems.step!(ds, simulation_time, true)

    lyap_scan[c,:] = [coupling_strength_range[c], DynamicalSystems.current_state(ds)[5]/DynamicalSystems.current_time(ds)]
end



axis_panel_d_inset = Axis(
    panel_d[1,1],
    width=Relative(0.55),
    height=Relative(0.55),
    halign=0.95,
    valign=0.95,
    limits = (-15,15, -0.5 , 0.1),
    xticks = [-20,-10,0,10,20],
    xtickformat = values -> [L"-20",L"-10",L"0",L"10",L"20"],
    xlabel = L"\alpha",
    ylabel = L"\lambda",
    xlabelpadding = -1,
    ylabelpadding = 0.5,
    xticklabelpad = 0.5,
    yticklabelpad = 1.5,
    yticks = [-0.4,-0.2, 0, 0.2],
    ytickformat = values -> [L"-0.4",L"-0.2",L"0",L"0.2"],
    xtickalign = 1,
    ytickalign = 1,    
    xgridvisible = false,
    ygridvisible = false
)

inset_plot_range = 1:(315*5)

lines!(
    axis_panel_d_inset,
    lyap_scan,
    color = default_color,
    linewidth = 1.5
)

lines!(
    axis_panel_d_inset,
    [-20.0,20.0],[0.0,0.0],
    color = TUD_grey80,
    linestyle = :dash,
    linewidth = 1
)

scatter!(
    axis_panel_d_inset,
    [-1],[-0.007010937880372775],
    color = default_color,
    markersize = 8
)


save(plotsdir("FIG3_limit_cycle_coupling_v2.pdf"), fig)
save(plotsdir("FIG3_limit_cycle_coupling_v2.png"), fig)


fig
