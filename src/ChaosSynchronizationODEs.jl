import StaticArrays 
import DifferentialEquations 

struct SystemParameters{Float64}
    couplingstrength::StaticArrays.SVector{3, Float64}
    numberofcouplingregions::Int64
    couplingpoints::Vector{StaticArrays.SVector{3, Float64}}
    couplingradii::Vector{Float64}
end





#ROESSLER

#a single Roessler attractor
function roesslerODE(x::StaticArrays.SVector{3, Float64}, p, t::Float64)
    dx1 = -x[2] - x[3]
    dx2 = x[1] + 0.2*x[2]
    dx3 = 0.2 + x[3]*(x[1] - 5.7)

    return StaticArrays.@SVector[dx1, dx2, dx3]
end

#two state-dependent master-slave-coupled Roessler attractors
function roesslerODE_coupled(x::StaticArrays.SVector{6, Float64}, p::SystemParameters, t::Float64) 
    #nonlinear master oscillator
    dx1 = - x[2] - x[3]
    dx2 = x[1] + 0.2*x[2]
    dx3 = 0.2 + x[3]*(x[1] - 5.7)

    #nonlinear slave-oscillator
    dx4 = - x[5] - x[6] 
    dx5 = x[4] + 0.2*x[5]
    dx6 = 0.2 + x[6]*(x[4] - 5.7)

    #coupling
    for k::Int in 1:p.numberofcouplingregions
        #if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 && (x[4]-p.couplingpoints[k][1])^2 + (x[5]-p.couplingpoints[k][2])^2 + (x[6]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= p.couplingstrength[1]*(x[4] - x[1])
                dx5 -= p.couplingstrength[2]*(x[5] - x[2])
                dx6 -= p.couplingstrength[3]*(x[6] - x[3])
                break
            end
        #end
    end

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6]
end

function roesslerODE_coupled_linear(x::StaticArrays.SVector{7, Float64}, p::SystemParameters, t::Float64)    
    #nonlinear master oscillator
    dx1 = - x[2] - x[3]
    dx2 = x[1] + 0.2*x[2]
    dx3 = 0.2 + x[3]*(x[1] - 5.7)

    #linearized slave-oscillator
    dx4 = - x[5] - x[6] 
    dx5 = x[4] + 0.2*x[5]
    dx6 = x[6]*(x[1] - 5.7) + x[3]*x[4]
    # dx7 is (v ⋅ dv/dt)
    dx7 = - x[4]*x[6] + 0.2*x[5]^2 + x[6]^2 * (x[1] - 5.7) + x[3]*x[4]*x[6]

    # coupling
    for k::Int in 1:p.numberofcouplingregions
        if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= p.couplingstrength[1]*x[4]
                dx5 -= p.couplingstrength[2]*x[5]
                dx6 -= p.couplingstrength[3]*x[6]
                dx7 -= p.couplingstrength[1]*x[4]^2 + p.couplingstrength[2]*x[5]^2 + p.couplingstrength[3]*x[6]^2
                break
            end
        end
    end
    
    # subtract radial component
    dx4 = dx4 - x[4] * dx7
    dx5 = dx5 - x[5] * dx7
    dx6 = dx6 - x[6] * dx7

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6, dx7]
end





#LORENZ

#a single Lorenz attractor
function lorenzODE(x::StaticArrays.SVector{3, Float64}, p, t::Float64)
    dx1 = 10*(x[2] - x[1])
    dx2 = x[1]*(28-x[3]) - x[2]
    dx3 = x[1]*x[2] - 8.0/3.0*x[3]

    return StaticArrays.@SVector[dx1, dx2, dx3]
end

#two state-dependent master-slave-coupled Lorenz attractors
function lorenzODE_coupled(x::StaticArrays.SVector{6, Float64}, p::SystemParameters, t::Float64) 
    #nonlinear master oscillator
    dx1 = 10*(x[2] - x[1])
    dx2 = x[1]*(28-x[3]) - x[2]
    dx3 = x[1]*x[2] - 8.0/3.0*x[3]

    #nonlinear slave-oscillator
    dx4 = 10*(x[5] - x[4])
    dx5 = x[4]*(28-x[6]) - x[5]
    dx6 = x[4]*x[5] - 8.0/3.0*x[6]

    #coupling
    for k::Int in 1:p.numberofcouplingregions
        #if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= p.couplingstrength[1]*(x[4] - x[1])
                dx5 -= p.couplingstrength[2]*(x[5] - x[2])
                dx6 -= p.couplingstrength[3]*(x[6] - x[3])
                break
            end
        #end
    end

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6]
end

function lorenzODE_coupled_linear(x::StaticArrays.SVector{7, Float64}, p::SystemParameters, t::Float64)    
    #nonlinear master oscillator
    dx1 = 10*(x[2] - x[1])
    dx2 = x[1]*(28-x[3]) - x[2]
    dx3 = x[1]*x[2] - 8.0/3.0*x[3]

    #linearized slave-oscillator
    dx4 = 10*(x[5] - x[4])
    dx5 = x[4]*(28-x[3]) - x[1]*x[6] - x[5]
    dx6 = x[4]*x[2] + x[1]*x[5] - 8.0/3.0*x[6]
    # dx7 is (v ⋅ dv/dt)
    dx7 = (10*(x[5] - x[4]))*x[4] + (x[4]*(28-x[3]) - x[5])*x[5] + (x[4]*x[2] - 8.0/3.0*x[6])*x[6]

    # coupling
    for k::Int in 1:p.numberofcouplingregions
        if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= p.couplingstrength[1]*x[4]
                dx5 -= p.couplingstrength[2]*x[5]
                dx6 -= p.couplingstrength[3]*x[6]
                dx7 -= p.couplingstrength[1]*x[4]^2 + p.couplingstrength[2]*x[5]^2 + p.couplingstrength[3]*x[6]^2
                break
            end
        end
    end
    
    # subtract radial component
    dx4 = dx4 - x[4] * dx7
    dx5 = dx5 - x[5] * dx7
    dx6 = dx6 - x[6] * dx7

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6, dx7]
end





#CHEN
# time rescaled to approximately match time scale of the dynamics witht the sampling time of other systems
# no qualitative changes, but all Lyapunov exponents must be scales down by the same factor of 10

#a single Chen attractor
function chenODE(x::StaticArrays.SVector{3, Float64}, p, t::Float64)
    dx1 = 0.1 * (35*(x[2] - x[1]))
    dx2 = 0.1 * (-7*x[1] - x[1]*x[3] + 28*x[2])
    dx3 = 0.1 * (x[1]*x[2] - 3*x[3])

    return StaticArrays.@SVector[dx1, dx2, dx3]
end

#two state-dependent master-slave-coupled Chen attractors
function chenODE_coupled(x::StaticArrays.SVector{6, Float64}, p::SystemParameters, t::Float64) 
    #nonlinear master oscillator
    dx1 = 0.1 * (35*(x[2] - x[1]))
    dx2 = 0.1 * (-7*x[1] - x[1]*x[3] + 28*x[2])
    dx3 = 0.1 * (x[1]*x[2] - 3*x[3])

    #nonlinear slave-oscillator
    dx4 = 0.1 * (35*(x[5] - x[4]))
    dx5 = 0.1 * (-7*x[4] - x[4]*x[6] + 28*x[5])
    dx6 = 0.1 * (x[4]*x[5] - 3*x[6])

    #coupling
    for k::Int in 1:p.numberofcouplingregions
        #if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= 0.1 * (p.couplingstrength[1]*(x[4] - x[1]))
                dx5 -= 0.1 * (p.couplingstrength[2]*(x[5] - x[2]))
                dx6 -= 0.1 * (p.couplingstrength[3]*(x[6] - x[3]))
                break
            end
        #end
    end

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6]
end

function chenODE_coupled_linear(x::StaticArrays.SVector{7, Float64}, p::SystemParameters, t::Float64)    
    #nonlinear master oscillator
    dx1 = 0.1 * (35*(x[2] - x[1]))
    dx2 = 0.1 * (-7*x[1] - x[1]*x[3] + 28*x[2])
    dx3 = 0.1 * (x[1]*x[2] - 3*x[3])

    #linearized slave-oscillator
    dx4 = 0.1 * (35*(x[5] - x[4]))
    dx5 = 0.1 * (-7*x[4] - x[1]*x[6] - x[4]*x[3] + 28*x[5])
    dx6 = 0.1 * (x[4]*x[2] + x[1]*x[5] - 3*x[6])
    # dx7 is (v ⋅ dv/dt)
    dx7 = 0.1 * ((35*(x[5] - x[4]))*x[4] + (-7*x[4] - x[3]*x[4] + 28*x[5])*x[5] + (x[4]*x[2] - 3*x[6])*x[6])

    # coupling
    for k::Int in 1:p.numberofcouplingregions
        if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= 0.1 * (p.couplingstrength[1]*x[4])
                dx5 -= 0.1 * (p.couplingstrength[2]*x[5])
                dx6 -= 0.1 * (p.couplingstrength[3]*x[6])
                dx7 -= 0.1 * (p.couplingstrength[1]*x[4]^2 + p.couplingstrength[2]*x[5]^2 + p.couplingstrength[3]*x[6]^2)
                break
            end
        end
    end
    
    # subtract radial component
    dx4 = dx4 - x[4] * dx7
    dx5 = dx5 - x[5] * dx7
    dx6 = dx6 - x[6] * dx7

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6, dx7]
end





#CHUA

#a single Chua attractor
function chuaODE(x::StaticArrays.SVector{3, Float64}, p, t::Float64)
    dx1 = 15.6 * (x[2]/5.0 - x[1] - (-5.0/7.0*x[1] - 3.0/14.0*(abs(x[1]+1)-abs(x[1]-1))))
    dx2 = 5.0 * 1.0 * (x[1] - x[2]/5.0 + x[3])
    dx3 = -28.0 * (x[2]/5.0)

    return StaticArrays.@SVector[dx1, dx2, dx3]
end

#two state-dependent master-slave-coupled Chua attractors
function chuaODE_coupled(x::StaticArrays.SVector{6, Float64}, p::SystemParameters, t::Float64) 
    #nonlinear master oscillator
    dx1 = 15.6 * (x[2]/5.0 - x[1] - (-5.0/7.0*x[1] - 3.0/14.0*(abs(x[1]+1)-abs(x[1]-1))))
    dx2 = 5.0 * 1.0 * (x[1] - x[2]/5.0 + x[3])
    dx3 = -28.0 * (x[2]/5.0)

    #nonlinear slave-oscillator
    dx4 = 15.6 * (x[5]/5.0 - x[4] - (-5.0/7.0*x[4] - 3.0/14.0*(abs(x[4]+1)-abs(x[4]-1))))
    dx5 = 5.0 * 1.0 * (x[4] - x[5]/5.0 + x[6])
    dx6 = -28.0 * (x[5]/5.0)

    #coupling
    for k::Int in 1:p.numberofcouplingregions
        #if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= p.couplingstrength[1]*(x[4] - x[1])
                dx5 -= p.couplingstrength[2]*(x[5] - x[2])
                dx6 -= p.couplingstrength[3]*(x[6] - x[3])
                break
            end
        #end
    end

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6]
end

function chuaODE_coupled_linear(x::StaticArrays.SVector{7, Float64}, p::SystemParameters, t::Float64)    
    #nonlinear master oscillator
    dx1 = 15.6 * (x[2]/5.0 - x[1] - (-5.0/7.0*x[1] - 3.0/14.0*(abs(x[1]+1)-abs(x[1]-1))))
    dx2 = 5.0 * 1.0 * (x[1] - x[2]/5.0 + x[3])
    dx3 = -28.0 * (x[2]/5.0)

    #linearized slave-oscillator
    dx4 = 15.6 * (x[5]/5.0 - x[4] - (-5.0/7.0*x[5] - 3.0/14.0*(sign(x[1]+1)*x[5]-sign(x[1]-1)*x[5])))
    dx5 = 5.0 * 1.0 * (x[4] - x[5]/5.0 + x[6])
    dx6 = -28.0 * (x[5]/5.0)
    # dx7 is (v ⋅ dv/dt)
    dx7 = (15.6 * (x[5]/5.0 - x[4] - (-5.0/7.0*x[5] - 3.0/14.0*(sign(x[1]+1)*x[5]-sign(x[1]-1)*x[5]))))*x[4] 
            + (5.0 * 1.0 * (x[4] - x[5]/5.0 + x[6]))*x[5] 
            + (-28.0 * (x[5]/5.0))*x[6]

    # coupling
    for k::Int in 1:p.numberofcouplingregions
        if( abs(x[1]-p.couplingpoints[k][1]) < p.couplingradii[k] && abs(x[2]-p.couplingpoints[k][2]) < p.couplingradii[k] && abs(x[3]-p.couplingpoints[k][3]) < p.couplingradii[k] )
            if( (x[1]-p.couplingpoints[k][1])^2 + (x[2]-p.couplingpoints[k][2])^2 + (x[3]-p.couplingpoints[k][3])^2 < p.couplingradii[k]^2 )
                dx4 -= p.couplingstrength[1]*x[4]
                dx5 -= p.couplingstrength[2]*x[5]
                dx6 -= p.couplingstrength[3]*x[6]
                dx7 -= p.couplingstrength[1]*x[4]^2 + p.couplingstrength[2]*x[5]^2 + p.couplingstrength[3]*x[6]^2
                break
            end
        end
    end
    
    # subtract radial component
    dx4 = dx4 - x[4] * dx7
    dx5 = dx5 - x[5] * dx7
    dx6 = dx6 - x[6] * dx7

    return StaticArrays.@SVector[dx1, dx2, dx3, dx4, dx5, dx6, dx7]
end





# -----------------------------------------------------------------------------
# NORMALIZATION CALLBACK
# -----------------------------------------------------------------------------

# Out-of-place (OOP) callback for immutable StaticArrays
function normalizecallback_affect!(integrator)
    u = integrator.u
    delta = sqrt(u[4]^2 + u[5]^2 + u[6]^2)
    v4 = u[4] / delta
    v5 = u[5] / delta
    v6 = u[6] / delta
    # Reconstruct the SVector
    integrator.u = StaticArrays.@SVector[u[1], u[2], u[3], v4, v5, v6, u[7]]
end
init_normalize = (c,u,t,integrator) -> normalizecallback_affect!(integrator) # Initialize state to be normalized
linearizedODE_normalizecallback = DifferentialEquations.DiscreteCallback((u,t,integrator) -> true, normalizecallback_affect!, initialize=init_normalize, save_positions=(false,false))

