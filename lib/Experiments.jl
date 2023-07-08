module Experiments

using ControlSystemsBase
using Plots, LaTeXStrings

using ControlSafetyBench
using ControlTimingSafety
using RealTimeScheduling

export nominal_trajectory_1, x_to_z_kill, x_to_z_skip_next, nominal_trajectory_2, slice_nominal, period_boosting, plotdev

function plotdev(sysc::AbstractStateSpace{<:Continuous}, 
        K::AbstractMatrix{<:Real}, p::Integer, newp::Integer, 
        z0::AbstractVecOrMat{<:Real}, safety_margin::Real, H::Integer, 
        az::Real, el::Real, ar::Real=1)
    # TODO: make nominal include full state, not just output (so we have 3-d plots)
    sysd = c2d(sysc, 0.001*p)
    
    a = hold_kill(sysd, K, MeetAny(1, 1))
    input = ones(Int64, H)
    nominal = evol(a, z0, input)'

    # Create empty plot with three dimensions
	plt = plot(xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t",
               camera=(az,el), aspect_ratio=ar)

    # Plot nominal trajectory
    plot!(nominal[1,:], nominal[2,:], 1:H+1, label="Nominal",
          color=:black, linewidth=1)

    # Plot safety pipe
	θ = LinRange(0, 2π, 40)
	circ_x = cos.(θ) * safety_margin
	circ_y = sin.(θ) * safety_margin
	for i in 1:H+1
		plot!(circ_x .+ nominal[1,i], circ_y .+ nominal[2,i], 
        repeat([i], size(θ,1)), label=(i == 1) ? "Safety Margin" : "", 
        seriestype=[:shape,], color=:lightblue)
	end

    plt
end

function period_boosting(sysc::AbstractStateSpace{<:Continuous}, 
        K::AbstractMatrix{<:Real}, newp::Integer, z0::AbstractVecOrMat{<:Real}, d_max::Real, maxwindow::Integer, n::Integer, nominal::Matrix{<:Real})
    sysd = c2d(sysc, 0.001*newp)
    nominal = slice_nominal(nominal, newp)
    H = size(nominal, 2) - 1

    synthesize_constraints(sysd, K, z0, d_max, maxwindow, n, H, fullresults=true, 
        nominal=nominal)[2]
end

function nominal_trajectory_1(sysc::AbstractStateSpace{<:Continuous},
        K::AbstractMatrix{<:Real}, x0::Real, p::Integer, H::Integer)
    sysd = c2d(sysc, 0.001*p)
    
    a = hold_kill(sysd, K, MeetAny(1, 1))
    z0 = x_to_z_kill(sysd, x0)
    input = ones(Int64, H)
    
    a.C * evol(a, z0, input)'
end

function nominal_trajectory_2(sysc::AbstractStateSpace{<:Continuous},
        K::AbstractMatrix{<:Real}, x0::Real, p::Integer, H::Integer)
    # Default with 1ms intervals
    sysd = c2d(sysc, 0.001)
    
    # No deadline requirements as we will be supplying the hit/miss pattern
    a = hold_skip_next(sysd, K)
    z0 = x_to_z_skip_next(sysd, x0)

    # Only the p-th deadline is hit
    input = Vector{Int64}((1:(H*p)) .% p .== 0)
    # Convert to 1=hit, 2=miss
    input = 2 .- input

    a.C * evol(a, z0, input)'
end

function x_to_z_kill(sys::AbstractStateSpace, x0::Real, u0::Real=0)
    vcat(repeat([Float64(x0)], sys.nx), repeat([Float64(u0)], sys.nu))
end

function x_to_z_skip_next(sys::AbstractStateSpace, x0::Real, u0::Real=0)
    vcat(repeat([Float64(x0)], sys.nx), repeat([0.0], sys.nx), repeat([Float64(u0)], sys.nu))
end

function slice_nominal(nom::AbstractMatrix{<:Real}, p::Integer)
    nom[:, (1:end) .% p .== 1]
end
    
end