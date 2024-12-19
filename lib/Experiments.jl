module Experiments

using ControlSystemsBase
using Plots
using LaTeXStrings

using ControlSafetyBench
using ControlTimingSafety
using RealTimeScheduling

export nominal_trajectory, trajectory, trajectory_coarse
export slice_nominal
export x_to_z_kill, x_to_z_skip_next
export period_boosting
export plotdev

function plotdev(nom::AbstractMatrix{<:Real}, nomp::Integer, x0::Real,
        d::Real, trj::Vector{<:Pair{String, <:AbstractMatrix{<:Real}}},
        az::Real, el::Real, ar::Real, legend::Symbol, title::String)
    t = size(nom, 2)
    @boundscheck all(x -> (x <= t), map(x -> size(x[2], 2), trj)) || 
        throw(ArgumentError("All trajectories must be of the same or shorter length "
        * "than the nominal trajectory."))
    
    limits = (-x0 - d*2, x0*1.2 + d*2)

    # Create empty plot with three dimensions
	plt = plot(xlabel=L"x position (meters)", ylabel=L"y position (meters)", zlabel=L"time (seconds)", legend=legend,
    title=title, camera=(az,el), aspect_ratio=ar, xlimits=limits, ylimits=limits)

    # Plot nominal trajectory
    plot!(nom[1,:], nom[2,:], (1:t)/nomp, label="Nominal", color=:black)

    # Plot safety pipe
	θ = LinRange(0, 2π, 20)
	circ_x = cos.(θ) * d
	circ_y = sin.(θ) * d
	for i in 1:nomp:t
		plot!(circ_x .+ nom[1,i], circ_y .+ nom[2,i], 
        repeat([i/nomp], size(θ,1)), label=(i == 1) ? "Safety Margin" : "", 
        seriestype=[:shape,], color=:lightblue)
	end

    for (label, γ) in trj
        plot!(γ[1,:], γ[2,:], (1:size(γ,2))/nomp, label=label)
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

nominal_trajectory(sysc::AbstractStateSpace{<:Continuous}, K::AbstractMatrix{<:Real},
        x0::Real, p::Integer, H::Integer) = trajectory(sysc, K, x0, p, p*H)

function trajectory(sysc::AbstractStateSpace{<:Continuous},
        K::AbstractMatrix{<:Real}, x0::Real, p::Integer, t::Integer;
        σ::Vector{<:Integer}=ones(Int64, t÷p))
    # Default with 1ms intervals
    sysd = c2d(sysc, 0.001)

    # No deadline requirements since this can see arbitrary input sequence.
    a = hold_skip_next(sysd, K)
    z0 = x_to_z_skip_next(sysd, x0)

    # Every p-th deadline is checked against the input sequence.
    # All other deadlines are misses
    σ = map(1:t) do x
        if x % p == 0 σ[Int64(x/p)] else 0 end
    end
    # Convert to 1=hit, 2=miss
    input = 2 .- σ

    # Since the states will be [x' x_saved' u']', we only want to extract x
    evol(a, z0, input)'[1:sysc.nx,:]
end

function trajectory_coarse(sysc::AbstractStateSpace{<:Continuous},
        K::AbstractMatrix{<:Real}, x0::Real, p::Integer, H::Integer;
        σ::Vector{<:Integer}=ones(Int64, H))
    sysd = c2d(sysc, 0.001*p)
    
    # No deadline requirements since this can see arbitrary input sequence.
    a = hold_kill(sysd, K)
    z0 = x_to_z_kill(sysd, x0)

    input = 2 .- σ

    # Since the states will be [x' u']', we only want to extract x
    evol(a, z0, input)'[1:sysc.nx,:]
end

"""
    x_to_z(sys, x0; skipnext=false)
    x_to_z(sys, x0, u0; skipnext=false)

Convert a given initial state `x0` (and optionally the initial input `u0`) to the 
initial state `z0` for the transducer automaton defined in `ControlTimingSafety`.
The `skipnext` flag determines if the state `z0` is in format for the `SkipNext` or
the `Kill` automaton.
"""
function x_to_z(sys::AbstractStateSpace, x0::Real, u0::Real=0; skipnext=false)
end

"""
    x_to_z_kill(sys, x0)
    x_to_z_kill(sys, x0, u0)

Convert a given initial state `x0` (with optional initial input `u0`) to
the required initial state `z0` for the `Kill` type Automaton.
"""
x_to_z_kill(sys::AbstractStateSpace, x0::Real, u0::Real=0) =
    @inbounds x_to_z_kill(sys, fill(x0, sys.nx), fill(u0, sys.nu))

x_to_z_kill(sys::AbstractStateSpace, x0::Vector{<:Real}, u0::Real=0) =
    x_to_z_kill(sys, x0, fill(u0, sys.nu))

x_to_z_kill(sys::AbstractStateSpace, x0::Real, u0::Vector{<:Real}=fill(0.0, sys.nu)) =
    x_to_z_kill(sys, fill(x0, sys.nx), u0)

function x_to_z_kill(sys::AbstractStateSpace, x0::Vector{<:Real}, u0::Vector{<:Real})
    @boundscheck length(x0) == sys.nx || 
        throw(ArgumentError("Dimension missmatch between `x0`` and `sys.nx`"))
    @boundscheck length(u0) == sys.nu || 
        throw(ArgumentError("Dimension missmatch between `u0`` and `sys.nu`"))
    Float64.(vcat(x0, u0))
end

"""
    x_to_z_skip_next(sys, x0)
    x_to_z_skip_next(sys, x0, u0)

Convert a given initial state `x0` (with optional initial input `u0`) to
the required initial state `z0` for the `SkipNext` type Automaton.
"""
x_to_z_skip_next(sys::AbstractStateSpace, x0::Real, u0::Real=0) =
    @inbounds x_to_z_skip_next(sys, fill(x0, sys.nx), fill(u0, sys.nu))

x_to_z_skip_next(sys::AbstractStateSpace, x0::Vector{<:Real}, u0::Real=0) =
    x_to_z_skip_next(sys, x0, fill(u0, sys.nu))

x_to_z_skip_next(sys::AbstractStateSpace, x0::Real, u0::Vector{<:Real}=fill(0.0, sys.nu)) =
    x_to_z_skip_next(sys, fill(x0, sys.nx), u0)

function x_to_z_skip_next(sys::AbstractStateSpace, x0::Vector{<:Real}, u0::Vector{<:Real})
    @boundscheck length(x0) == sys.nx || 
        throw(ArgumentError("Dimension missmatch between `x0`` and `sys.nx`"))
    @boundscheck length(u0) == sys.nu || 
        throw(ArgumentError("Dimension missmatch between `u0`` and `sys.nu`"))
    Float64.(vcat(x0, fill(0.0, sys.nx), u0))
end

function slice_nominal(nom::AbstractMatrix{<:Real}, p::Integer)
    nom[:, (1:end) .% p .== 1]
end
    
end