using LinearAlgebra
using ControlSystemsBase
using Distances
using QuadGK
using JLD
using DelimitedFiles

# using SparseArrays
# using OffsetArrays
# using Distributions

struct Automaton
	# L: locations.  Legal locations are in the range 1:L.
	L::Int64
	# A: actions (input alphabet).  Legal actions are in the range 1:A.
	A::Int64
	# Φ: dynamics matrices (output alphabet).  Array of square matrices of equal size.
	Φ::Vector{AbstractMatrix{Float64}}
	# T: transition function.  T[l,a] is a location in 1:L, or missing.
	T::Matrix{Union{Missing, Int64}}
	# μ: output function.  μ[l,a] is a character from Φ.
	μ::Matrix{Union{Missing, Int64}}
	# l_int: initial location in L.
	l_int::Int64
end

function Automaton_lint(other::Automaton, l_int::Int64)
    Automaton(other.L, other.A, other.Φ, other.T, other.μ, l_int)
end

function Evol(z_0, automaton, input)
    t_max = size(input, 1)
    z = zeros(size(z_0, 1), t_max + 1)
    z[:,1] = z_0
	l = automaton.l_int
	# For each time step
    for t = 1:t_max
		# Get the dynamics matrix
		μ = automaton.μ[l, input[t]]
		# If we hit a missing transition, return the states that we reached,
		# and a missing final location to signal the problem to the caller.
		if ismissing(μ)
			return z[:,1:t]', missing
		end
		# Apply the dynamics
        #z[:,t+1] = automaton.Φ[:,:,μ] * z[:,t]
		z[:,t+1] = automaton.Φ[μ] * z[:,t]
		# Transition to the new location
		l = automaton.T[l, input[t]]
    end
    z', l
end

function Augment(x, automaton)
    [x; zeros(size(automaton.Φ[1],1) - size(x,1), size(x,2))]
end

function HoldAndKill(sysd, K, miss=nothing, window=1)
	p, r = size(sysd.B)
	
	# Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end
	
	# Define the automaton's structure
	if miss === nothing || miss == -1
		L = 1
		T = [1 1]
		μ = [2 1]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		min_hits = window - miss
		L = 2^window
		T = zeros(Union{Missing, Int64}, (L, 2))
		μ = zeros(Union{Missing, Int64}, (L, 2))
		for i = 0:L-1
			T[i+1, 1] = ((i << 1) & (L - 1)) + 1		# miss
			T[i+1, 2] = ((i << 1) & (L - 1) | 1) + 1	# hit
			μ[i+1, 1] = 2
			μ[i+1, 2] = 1
			for j = 1:2
				if count_ones(T[i+1, j] - 1) < min_hits
					T[i+1, j] = missing
					μ[i+1, j] = missing
				end
			end
		end
		T = L + 1 .- T
		reverse!(T, dims=1)
		reverse!(μ, dims=1)
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
		    # Φ_H
			[sysd.A  sysd.B;
			 K_x  K_u],
			# Φ_M
			[sysd.A  sysd.B;
	         zeros(r, p)  I(r)]],
		T, μ, 1)
end

function ZeroAndKill(sysd, K, miss=nothing, window=1)
	p, r = size(sysd.B)
	
	# Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end
	
	# Define the automaton's structure
	if miss === nothing || miss == -1
		L = 1
		T = [1 1]
		μ = [2 1]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		min_hits = window - miss
		L = 2^window
		T = zeros(Union{Missing, Int64}, (L, 2))
		μ = zeros(Union{Missing, Int64}, (L, 2))
		for i = 0:L-1
			T[i+1, 1] = ((i << 1) & (L - 1)) + 1		# miss
			T[i+1, 2] = ((i << 1) & (L - 1) | 1) + 1	# hit
			μ[i+1, 1] = 2
			μ[i+1, 2] = 1
			for j = 1:2
				if count_ones(T[i+1, j] - 1) < min_hits
					T[i+1, j] = missing
					μ[i+1, j] = missing
				end
			end
		end
		T = L + 1 .- T
		reverse!(T, dims=1)
		reverse!(μ, dims=1)
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
		    # Φ_H
			[sysd.A  sysd.B;
			 K_x  K_u],
			# Φ_M
			[sysd.A  sysd.B;
	         zeros(r, p + r)]],
		T, μ, 1)
end

function HoldAndSkipNext(sysd, K, miss=nothing, window=1)
	p, r = size(sysd.B)
	
	# Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end
	
	# Define the automaton's structure
	if miss === nothing || miss == -1
		# No constraint. Any # of misses are allowed.
		L = 2
		T = [2 1;
			 2 1]
		μ = [3 1;
			 4 2]
	elseif miss == 0
		# No miss allowed
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		min_hits = window - miss
		L = 2^window
		T = zeros(Union{Missing, Int64}, (L, 2))
		μ = zeros(Union{Missing, Int64}, (L, 2))
		for i = 0:L-1
			T[i+1, 1] = ((i << 1) & (L - 1)) + 1		# miss
			T[i+1, 2] = ((i << 1) & (L - 1) | 1) + 1	# hit
			μ[i+1, 1] = 4 - (i & 1)
			μ[i+1, 2] = 2 - (i & 1)
			for j = 1:2
				if count_ones(T[i+1, j] - 1) < min_hits
					T[i+1, j] = missing
					μ[i+1, j] = missing
				end
			end
		end
		T = L + 1 .- T
		reverse!(T, dims=1)
		reverse!(μ, dims=1)
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
			# Φ_HH
			[sysd.A  zeros(p, p)  sysd.B;
			 zeros(p, 2p + r);
			 K_x  zeros(r, p)  K_u],
			# Φ_MH
			[sysd.A  zeros(p, p)  sysd.B;
	         zeros(p, 2p + r);
	         zeros(r, p)  K_x  K_u],
			# Φ_HM
			[sysd.A  zeros(p, p)  sysd.B;
	         I(p)  zeros(p, p + r);
	         zeros(r, 2p)  I(r)],
			# Φ_MM
			[sysd.A  zeros(p, p)  sysd.B;
	         zeros(p, p)  I(p)  zeros(p, r);
	         zeros(r, 2p)  I(r)]],
		T, μ, 1)
end

function ZeroAndSkipNext(sysd, K, miss=nothing, window=1)
	p, r = size(sysd.B)
	
	# Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end
	
	# Define the automaton's structure
	if miss === nothing || miss == -1
		# No constraint. Any # of misses are allowed.
		L = 2
		T = [2 1;
			 2 1]
		μ = [3 1;
			 4 2]
	elseif miss == 0
		# No misses alowed
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		min_hits = window - miss
		L = 2^window
		T = zeros(Union{Missing, Int64}, (L, 2))
		μ = zeros(Union{Missing, Int64}, (L, 2))
		for i = 0:L-1
			T[i+1, 1] = ((i << 1) & (L - 1)) + 1		# miss
			T[i+1, 2] = ((i << 1) & (L - 1) | 1) + 1	# hit
			μ[i+1, 1] = 4 - (i & 1)
			μ[i+1, 2] = 2 - (i & 1)
			for j = 1:2
				if count_ones(T[i+1, j] - 1) < min_hits
					T[i+1, j] = missing
					μ[i+1, j] = missing
				end
			end
		end
		T = L + 1 .- T
		reverse!(T, dims=1)
		reverse!(μ, dims=1)
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
		    # Φ_HH
			[sysd.A  zeros(p, p)  sysd.B;
			 zeros(p, 2p + r);
			 K_x  zeros(r, p)  K_u],
			# Φ_MH
			[sysd.A  zeros(p, p)  sysd.B;
	         zeros(p, 2p + r);
	         zeros(r, p)  K_x  K_u],
			# Φ_HM
			[sysd.A  zeros(p, p)  sysd.B;
	         I(p)  zeros(p, p + r);
	         zeros(r, 2p + r)],
			# Φ_MM
			[sysd.A  zeros(p, p)  sysd.B;
	         zeros(p, p)  I(p)  zeros(p, r);
	         zeros(r, 2p + r)]],
		T, μ, 1)
end

strat_map = Dict(
    "HK" => HoldAndKill,
    "ZK" => ZeroAndKill,
    "HS" => HoldAndSkipNext,
    "ZS" => ZeroAndSkipNext
)
strat_names = sort([keys(strat_map)...])

function c2da(sysc::AbstractStateSpace{<:Continuous}, Ts::Real, d::Real)
    f_Φ = s -> ℯ^(sysc.A * s)
    Φ = f_Φ(Ts)
    Γ_0 = quadgk(f_Φ, 0, Ts-d)[1] * sysc.B
    Γ_1 = quadgk(f_Φ, Ts-d, Ts)[1] * sysc.B
    Φ_a = [Φ  Γ_1;
		   zeros(size(sysc.B, 2), size(sysc.A, 2)+size(sysc.B, 2))]
    Γ_a = [Γ_0;
		   I]
    C_a = [sysc.C  zeros(size(sysc.C,1), size(sysc.B,2))]
    D_a = zeros(size(C_a,1), size(Γ_a,2))
    ss(Φ_a, Γ_a, C_a, D_a, Ts)
end

function BoundedTreeFast(automaton, bounds, n)
	# Stack
	z = Matrix{Float64}(undef, n+1, size(bounds, 1))
	loc = Vector{Int64}(undef, n+1)
	act = Vector{Int64}(undef, n+1)
	
	# Bounding boxes for each final location, time step
	ret = Array{Float64}(undef, automaton.L, n+1, size(bounds, 1), 2) * NaN

	# For each corner
	for corner in Base.product(eachrow(bounds)...)
		# Create the stack frame for time 0
		z[1,:] = [c for c in corner]
		loc[1] = automaton.l_int
		act[1] = 1
		# Initialize the stack pointer
		sp = 1
		# While we haven't popped all the way out
		while sp >= 1
			# If we've reached a leaf
			if sp == n+1
				# Calculate min and max for this final location at each time step
				ret[loc[sp],:,:,1] = minimum(x->(isnan(x) || isinf(x)) ? Inf : x, cat(z, ret[loc[sp],:,:,1], dims=3), dims=3)
				ret[loc[sp],:,:,2] = maximum(x->(isnan(x) || isinf(x)) ? -Inf : x, cat(z, ret[loc[sp],:,:,2], dims=3), dims=3)
				sp -= 1
			# If we're out of actions from this step
			elseif act[sp] > automaton.A
				sp -= 1
			# If the transition is missing
			elseif ismissing(automaton.T[loc[sp], act[sp]])
				# Try the next transition
				act[sp] = act[sp] + 1
			# If the transition is present
			else
				z[sp+1,:] = automaton.Φ[automaton.μ[loc[sp], act[sp]]] * z[sp,:]
				loc[sp+1] = automaton.T[loc[sp], act[sp]]
				act[sp+1] = 1
				act[sp] = act[sp] + 1
				sp = sp + 1
			end
		end
	end
	ret
end

function get_traj(automaton, bounds, steps; dims=axes(bounds,1), hitpattern=repeat([2], steps-1))
	# Dimensions: state variables, points, time
	ev = Array{Float64}(undef, length(dims), 2^size(bounds,1), steps)
	corners = corners_from_bounds(bounds, dims=axes(bounds,1))
	for (i, c) in enumerate(eachcol(corners))
		e, _ = Evol(c, automaton, hitpattern)
		ev[:,i,:] = e'[dims,:]

		# For DATE 23, we currently only consider one initial location (instead of an interval)
		break
		# ========
	end
	
	# For DATE 23, we currently only consider one initial location (instead of an interval)
	ev[:,1,:]'
	# ========
end

function deviation(automaton, bounds, reachable; dims=axes(bounds,1), metric=Euclidean(), nominal=repeat([2],size(reachable,1)-1))
	if nominal === nothing
		# Assume a nominal behavior of all hits
		nominal = ones(Int64, size(reachable, 1) - 1) .+ 1
	end
	
	if dims === nothing
		# Assume dimensions 1 and 2 are the plant state
		# TODO: this is mostly just hardcoded for now, ignoring the dims variable.
		#   Partially the fault of corners_from_bounds.
		dims = (1, 2)
	end
	
	# Dimensions: state variables, points, time
	reachable_corners = cat([corners_from_bounds(reachable[t,:,:], dims=dims) for t in axes(reachable, 1)]..., dims=3)
	
	# Dimensions: state variables, points, time
	ev = Array{Float64}(undef, length(dims), 2^size(bounds,1), size(reachable, 1))
	corners = corners_from_bounds(bounds, dims=axes(bounds,1))
	for (i, c) in enumerate(eachcol(corners))
		e, _ = Evol(c, automaton, nominal)
		ev[:,i,:] = e'[dims,:]
	end
	
	# Compute Hausdorff distance at each time step
	H = Array{Float64}(undef, size(ev, 3))
	for t in axes(ev, 3)
		dist = pairwise(metric, reachable_corners[:,:,t], ev[:,:,t])
		H_row = maximum(minimum(dist, dims=1))
		H_col = maximum(minimum(dist, dims=2))
		H[t] = maximum((H_row, H_col))
	end
	H
end

function BoundedTreeIter(automaton, bounds, n, t, safety_margin; dims=axes(bounds,1))::Tuple{Float64, Int64}
	q = size(bounds, 1)
	
	# Dimensions: time, augmented state, min/max
	all_bounds = Array{Float64}(undef, n*t+1, q, 2)
	all_bounds[1,:,:] = bounds
	
	bounds = BoundedTreeFast(automaton, bounds, n)
	all_bounds[1:n+1,:,:] = merge_bounds(bounds)[:,:,1:end]

	temp = deviation(automaton, all_bounds[1,:,:], all_bounds[1:n+1,:,:], dims=dims)
	if maximum(temp) > safety_margin
		return (maximum(temp), argmax(temp))
	end
	
	# Dimensions: initial location, final location, time, augmented state, min/max
	new_bounds = Array{Any}(undef, automaton.L, automaton.L, n+1, q, 2)
	for i in 1:t-1
		# Simulate each box from previous iteration
		for i in 1:automaton.L
			a = Automaton_lint(automaton, i)
			new_bounds[i,:,:,:,:] = BoundedTreeFast(a, bounds[i,end,:,:], n)
		end
		# Merge resulting boxes from these simulations
		for i in 1:automaton.L
			bounds[i,:,:,:] = merge_bounds(new_bounds[:,i,:,:,:])
		end
		# Save the bounds
		all_bounds[n*i+2:n*(i+1)+1,:,:] = merge_bounds(bounds)[2:end,:,:]
		
		temp = deviation(automaton, all_bounds[1,:,:], all_bounds[1:n*(i+1)+1,:,:], dims=dims)
		if maximum(temp) > safety_margin
			return (maximum(temp), argmax(temp))
		end
	end

	(maximum(temp), argmax(temp))
end

function merge_bounds(b)
	mins = minimum(x->(isnan(x) || isinf(x)) ? Inf : x, b[:,:,:,1], dims=1)
	maxs = maximum(x->(isnan(x) || isinf(x)) ? -Inf : x, b[:,:,:,2], dims=1)
	cat(dims=4, mins, maxs)[1,:,:,:]
end

function corners_from_bounds(bounds; cycle=false, dims=nothing)
	if dims === nothing
		dims = axes(bounds, 1)
	end
	ldims = length(dims)

	corners = cat(reshape([[c...] for c in Base.product(eachrow(bounds[dims,:])...)], 2^ldims)..., dims=2)
	if cycle
		gray(x) = x ⊻ (x >> 1)
		[corners[:,gray.(0:2^ldims-1) .+ 1] corners[:,1]]
	else
		corners
	end
end