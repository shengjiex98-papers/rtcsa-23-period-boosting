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

function plotdev(sysc::AbstractStateSpace{<:Continuous}, 
        nom::AbstractMatrix{<:Real}, nomp::Integer, safety_margin::Real,
        trj::Vector{Pair{String, <:AbstractMatrix{<:Real}}};
        az::Real=50, el::Real=35, ar::Real=1)
    t = size(nom, 2)
    @boundscheck all(isequal(t), map(x -> size(x[2], 2), trj)) || 
        throw(ArgumentError("All trajectories must be the same length."))

    # Create empty plot with three dimensions
	plt = plot(xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t",
               camera=(az,el), aspect_ratio=ar)

    # Plot nominal trajectory
    plot!(nominal[1,:], nominal[2,:], (1:t)/1000, label="Nominal", color=:black)

    # Plot safety pipe
	θ = LinRange(0, 2π, 20)
	circ_x = cos.(θ) * safety_margin
	circ_y = sin.(θ) * safety_margin
	for i in 1:nomp:safety_margin
		plot!(circ_x .+ nominal[1,i], circ_y .+ nominal[2,i], 
        repeat([i], size(θ,1)), label=(i == 1) ? "Safety Margin" : "", 
        seriestype=[:shape,], color=:lightblue)
	end

    for (label, γ) in trj
        plot(γ[1,:], γ[2,:], (1:t)/1000, label=label)
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