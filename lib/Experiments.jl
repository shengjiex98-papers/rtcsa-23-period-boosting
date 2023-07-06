module Experiments

using ControlSystemsBase

using ControlSafetyBench
using ControlTimingSafety
using RealTimeScheduling

export nominal_trajectory_1, x_to_z_kill, x_to_z_skip_next, nominal_trajectory_2, slice_nominal, period_boosting

function period_boosting(sysc::AbstractStateSpace{<:Continuous}, K::AbstractMatrix{<:Real},
        newp::Integer, z0::AbstractVecOrMat{<:Real}, d_max::Real, maxwindow::Integer, 
        n::Integer, nominal::Matrix{<:Real})
    sysd = c2d(sysc, 0.001*newp)
    nominal = slice_nominal(nominal, newp)
    H = size(nominal, 1) - 1

    synthesize_constraints(sysd, K, z0, d_max, maxwindow, n, H, fullresults=true, 
        nominal=nominal)[2]
end

function nominal_trajectory_1(sysc::AbstractStateSpace{<:Continuous},
        x0::Real, p_ms::Integer, H::Integer)
    sysd = c2d(sysc, 0.001*p_ms)
    K = delay_lqr(sysc, 0.001*p_ms)
    
    a = hold_kill(sysd, K, MeetAny(1, 1))
    z0 = x_to_z_kill(sysd, x0)
    input = ones(Int64, H)
    
    evol(a, z0, input)
end

function nominal_trajectory_2(sysc::AbstractStateSpace{<:Continuous},
        x0::Real, p_ms::Integer, H::Integer)
    # Default with 1ms intervals
    sysd = c2d(sysc, 0.001)
    K = delay_lqr(sysc, 0.001*p_ms)
    
    # No deadline requirements as we will be supplying the hit/miss pattern
    a = hold_skip_next(sysd, K)
    z0 = x_to_z_skip_next(sysd, x0)

    # Only the p_ms-th deadline is hit
    input = Vector{Int64}((1:(H*p_ms)) .% p_ms .== 0)
    # Convert to 1=hit, 2=miss
    input = 2 .- input

    evol(a, z0, input)[:, vcat(1:sysc.nx, 2*sysc.nx+1:2*sysc.nx+sysc.nu)]
end

function x_to_z_kill(sys::AbstractStateSpace, x0::Real, u0::Real=0)
    vcat(repeat([Float64(x0)], sys.nx), repeat([Float64(u0)], sys.nu))
end

function x_to_z_skip_next(sys::AbstractStateSpace, x0::Real, u0::Real=0)
    vcat(repeat([Float64(x0)], sys.nx), repeat([0.0], sys.nx), repeat([Float64(u0)], sys.nu))
end

function slice_nominal(nom::AbstractMatrix{<:Real}, p_ms::Integer)
    nom[(1:end) .% p_ms .== 1,:]
end
    
end