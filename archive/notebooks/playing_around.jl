### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fd140d06-742c-44a9-87ee-304533065222
begin
	using Pkg
	
	using LinearAlgebra
	using ControlSystems
	using Random
	using Plots
	gr()
	using PlutoUI
	using LaTeXStrings
	using QuadGK
	using Distances
	using SparseArrays
	using OffsetArrays
	using Distributions

	using CSV
	using Tables
	using Printf
	using Markdown

	md"""Import needed packages."""
end

# ╔═╡ ee503656-cd73-11ec-13cb-61d08f418ed6
md"""
# Synthesizing Weakly-Hard Bounds

We are given a dynamical system and a trajectory in its
state space resulting from its dynamics. We refer to this as the
nominal trajectory, construct a safety pipe around it, and label
all trajectories in this pipe as safe. All such safe trajectories
constitute a property 𝜙. 

When the controller of this dynamical
system is implemented as software running on a computational
platform, it suffers from deadline hits and misses because of
other tasks sharing the platform. Such deadline hits and misses
result in the system deviating from its nominal trajectory. The
question we are interested in asking is:
Which deadline hit/miss patterns result in all
trajectories satisfying 𝜙?

## Automata basics

To begin our implementation, we first define an automaton struct.
"""

# ╔═╡ a402d61b-d832-42b1-83f8-de95a87cd771
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

# ╔═╡ 39568986-8cf7-4fe2-81b9-8139f26f14fa
md"""
For convenience later on, we'll also want this little function to make a new Automaton that's an exact copy of another, but with a different initial location.
"""

# ╔═╡ 110ca59b-f6e1-4033-9570-8e2b135f4f93
Automaton_lint(other::Automaton, l_int::Int64) = Automaton(other.L, other.A, other.Φ, other.T, other.μ, l_int)

# ╔═╡ 3dc10fd0-9a52-4090-a3d5-4ff048c5a28b
md"""
Now we define a function that runs an automaton, given an initial augmented state and string of actions.  It returns a matrix whose columns are the state vector at each time step, and the final location in the automaton.  This information allows us to continue the evolution in a second call to Evol, which we will use later.
"""

# ╔═╡ edd5649e-efcc-49bf-a333-0e14e8be518e
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

# ╔═╡ d1f0cdfb-f6b6-4bac-a3e1-864be6e9be5e
md"""
One more utility function next, to pad a state $x$ (or array of states) into an augmented state $z$ used by an automaton.  This will be convenient when we use the automata later on.
"""

# ╔═╡ e2a404c9-ed92-4367-bad2-f425606b88dc
Augment(x, automaton) = [x; zeros(size(automaton.Φ[1],1) - size(x,1), size(x,2))]

# ╔═╡ 8025bd64-5956-455b-b232-e41cb541ca09
md"""
We implement a function to construct a Hold&Skip-Next automaton, given a discrete state-space model and a feedback gain matrix K. Under this strategy, on a deadline miss, the old control input is held, while the job that missed is allowed to keep running and the next jobs are skipped until it completes.

!!! note
	The semantics of Hold&Skip-Next are following [^maggiostability], but the implementation is different, allowing for any number of sequential deadline misses.  This makes the augmented matrices much smaller, resulting in much faster computation.
"""

# ╔═╡ 8e910c95-76dc-45c1-b09f-b0b725e49b7d
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

# ╔═╡ ec833777-4b93-48e7-a076-dbcf147a908a
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

# ╔═╡ 41162c87-f54f-4e8c-8cfd-ac7a3916ec3d
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

# ╔═╡ 0216d280-02e6-4c24-b6e6-e1ff7b716e08
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

# ╔═╡ fb97808c-0911-4d19-a29d-9ab97e9d3c5f
md"""
## System models

Next, we present a few dynamical systems for our study. They are: RC nework, electrical steering, aircraft pitch, and F1tenth car.
"""

# ╔═╡ fd65777a-0e7e-41a7-a11a-f6545bf9cef4
begin
	sys_rc = let
		r_1 = 100000
		r_2 = 500000
		r_3 = 200000
		c_1 = 0.000002
		c_2 = 0.000010
		A = [-1/c_1 * (1/r_1 + 1/r_2)  1/(r_2*c_1)
		     1/(r_2*c_2)               -1/c_2 * (1/r_2 + 1/r_3)]
		B = [1/(r_1*c_1)
		     1/(r_3*c_2)]
		C = [1  -1]
		D = 0
		ss(A, B, C, D)
	end
	sys_es = let
		R = 0.025
		ωel = 2000*pi
		Ld = 0.0001
		Lq = 0.00012
		A = [-R/Ld  Lq * ωel / Ld;  -Ld * ωel / Lq   -R / Lq]
		B = [1 / Ld   0;  0   1/Lq]
		C = [0 0]
		D = [0 0]
		
		ss(A, B, C, D)
	end
	sys_ap = let
		A = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
		B = [0.232; 0.0203; 0];
		C = [0 0 1];
		D = [0];
	
		ss(A, B, C, D)
	end
	sys_rc
end

# ╔═╡ 083caada-cb9b-422f-9f1d-e72ce127f66e
md"""
And the discrete versions of each model.
"""

# ╔═╡ a6babf6a-48b0-4ff6-9b43-7b9871f506c8
sysd_rc = c2d(sys_rc, 0.02)

# ╔═╡ ce6cd63a-6393-4a3b-92c0-6ed3fb5edd79
sysd_es = c2d(sys_es, 0.02)

# ╔═╡ 713a1e90-052d-42a7-b61d-8001da3c2129
sysd_ap = c2d(sys_ap, 0.02)

# ╔═╡ 83347d7d-b559-4ee1-9272-2a249a378113
sysd_f1 = let
	A = [1 0.13; 0 1]
	B = [0.02559; 0.3937]
	C = [0 0]
	D = 0
	ss(A, B, C, D, 0.02)
end

# ╔═╡ de96915a-decc-4193-a63a-9964707dc025
sysd_us = let
	A = [1.1053 0; -0.0209 0.99]
	B = [0.0526 0.0105; -0.0393 0.0994]
	C = [1 0]
	D = [0 0]
	ss(A, B, C, D, 0.1)
end

# ╔═╡ dc7b7a0e-3b0e-4c32-bc8e-9b1fe319702c
sysd_me = let
	A = [0 1; -0.8 1.8]
	B = [0; 1]
	C = [0 1]
	D = 0
	ss(A, B, C, D, 0.1)
end

# ╔═╡ 839c0a6c-f6c2-4f53-abae-63739b4ad75b
md"""
Finally we introduce delay of one period to each model and compute a controller K.
"""

# ╔═╡ 9ff9233d-d0b9-4687-a0d5-4a396a8bad54
function c2da(sysc::AbstractStateSpace{<:ControlSystems.Continuous}, Ts::Real, d::Real)
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

# ╔═╡ ec282dbc-e4e3-4197-8587-130b34ee9550
md"""
We now create an optimal controller for the delayed system using LQR.
"""

# ╔═╡ 805e54ef-be8a-4956-80a0-700ec1b96a92
K_rc = let
	Q = I * 2
	R = I * 1
	sysd_rc_delay = c2da(sys_rc, 0.02, 0.02)
	lqr(sysd_rc_delay, Q, R)
end

# ╔═╡ 01dacf35-a7d4-4e2b-8391-d40fe3b317f7
K_es = let
	Q = I
	R = I
	[lqr(sysd_es, Q, R) zeros(2, 2)]
end

# ╔═╡ 974d1e72-c6e8-4d61-acc8-a5263307bd51
K_ap = let
	Q = [0 0 0 0; 0 0 0 0; 0 0 50 0; 0 0 0 0]
	R = I
	sysd_ap_delay = c2da(sys_ap, 0.02, 0.02)
	lqr(sysd_ap_delay, Q, R)
end

# ╔═╡ a43c4d8e-1935-4225-8db0-139b7b4207a6
K_f1 = [0.293511 0.440267 0]

# ╔═╡ 96a0379f-ffb3-47b4-8719-d60ada971d42
K_us = [4.7393 0.2430; 0.2277 -0.8620]

# ╔═╡ b285273a-04fd-494f-bb15-ab0463e07496
K_me = [0.256 -0.372 0]

# ╔═╡ f6ebc370-9059-4236-836e-d080bd012fb9
begin
	model_map = Dict(
		"RC" => (sysd_rc, K_rc),
		"ES" => (sysd_es, K_es),
		"AP" => (sysd_ap, K_ap),
		"F1" => (sysd_f1, K_f1),
		"US" => (sysd_us, K_us),
		"ME" => (sysd_me, K_me)
	)
	model_names = sort([keys(model_map)...])
end

# ╔═╡ 731f357a-a4c5-44dd-85c2-105a950d72e6
md"""
## Bounded Tree Reachable Set

In this section, we use a *bounded tree method* to compute a reachable set for the system under any sequence of deadline hits and misses.  In this method, we start with an axis-aligned box initial set, and simulate its evolution for $n$ steps, for each of the $|\mathcal{A}|^n$ possible sequences of scheduler actions.  We then take an axis-aligned bounding box of the resulting sets, and repeat this process.

This algorithm more generally supports running any automaton, including those that implement different deadline miss strategies, and weakly-hard constraints.  This gives us more power to model scheduler behavior.

We begin by defining a function to run one iteration of this process.  For the iteration to be sound, each iteration must return one bounding box for each final location in the automaton.  Then, a separate tree will be run for each final location output by the previous round, so this function must also take as a parameter the initial location to use. 

This algorithm leverages a preorder traversal of the automaton to remove the redundant computation, bringing the time complexity down to $O(|\mathcal{A}|^n)$.
"""

# ╔═╡ 811d1e05-a6ca-43bd-9f62-395392004747
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

# ╔═╡ 48bfd3c8-87ba-482d-b0e2-22019e1afef4
md"""
Here we need to define two helper functions.
"""

# ╔═╡ c2edd2a4-0ad1-45dc-8d6e-4c922a927e3c
function merge_bounds(b)
	mins = minimum(x->(isnan(x) || isinf(x)) ? Inf : x, b[:,:,:,1], dims=1)
	maxs = maximum(x->(isnan(x) || isinf(x)) ? -Inf : x, b[:,:,:,2], dims=1)
	cat(dims=4, mins, maxs)[1,:,:,:]
end

# ╔═╡ e014f895-6ed1-4489-8288-68e184f2ba39
function BoundedTreeIterOld(automaton, bounds, n, t)
	q = size(bounds, 1)
	
	# Dimensions: time, augmented state, min/max
	all_bounds = Array{Float64}(undef, n*(t+1)+1, q, 2)
	all_bounds[1,:,:] = bounds
	
	bounds = BoundedTreeFast(automaton, bounds, n)
	all_bounds[1:n+1,:,:] = merge_bounds(bounds)[:,:,1:end]
	
	# Dimensions: initial location, final location, time, augmented state, min/max
	new_bounds = Array{Any}(undef, automaton.L, automaton.L, n+1, q, 2)
	for i in 1:t
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
	end
	all_bounds
end

# ╔═╡ 515802ce-3186-43d9-9897-6286622782c2
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

# ╔═╡ a32a9411-1e3b-47dc-81f7-89ccf851d683
md"""
## Computing an Upper-Bound on Deviation

Given reachable sets at each time step, it is straightforward to compute a bound to the deviation from a nominal trajectory.  We can simply compute the system evolution from the corners of the initial set using the nominal behavior, then calculate the distance between these corners at each time step and the corners of the corresponding reachable set.  Taking the maximum yields our deviation bound $\mathbf{d}_\mathit{ub}$.
"""

# ╔═╡ 550655ed-411c-4fbc-b932-730bc32514fc
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

# ╔═╡ 8c99c664-3cea-40a0-83c2-d8d1f6c8e46f
function BoundedTreeIter(automaton, bounds, n, t, deviation_bound; dims=axes(bounds,1))::Tuple{Float64, Int64}
	q = size(bounds, 1)
	
	# Dimensions: time, augmented state, min/max
	all_bounds = Array{Float64}(undef, n*t+1, q, 2)
	all_bounds[1,:,:] = bounds
	
	bounds = BoundedTreeFast(automaton, bounds, n)
	all_bounds[1:n+1,:,:] = merge_bounds(bounds)[:,:,1:end]

	temp = deviation(automaton, all_bounds[1,:,:], all_bounds[1:n+1,:,:], dims=dims)
	if maximum(temp) > deviation_bound
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
		if maximum(temp) > deviation_bound
			return (maximum(temp), argmax(temp))
		end
	end

	(maximum(temp), argmax(temp))
end

# ╔═╡ 3873d023-15f1-490c-b3b7-caf6b96f8301
# ╠═╡ disabled = true
#=╠═╡
function violate(dims, reachable, ev, deviation_bound; metric=Euclidean())
	reachable_corners = cat([corners_from_bounds(reachable[t,:,:], dims=dims) for t in axes(reachable, 1)]..., dims=3)
	for t in axes(ev, 3)
		dist = pairwise(metric, reachable_corners[:,:,t], ev[:,:,t])
		H_row = maximum(minimum(dist, dims=1))
		H_col = maximum(minimum(dist, dims=2))
		H = maximum((H_row, H_col))
		if H > deviation_bound
			return true
		end
	end
	return false
end
  ╠═╡ =#

# ╔═╡ 80174b73-34f9-4c62-bf09-29672b68a558
md"""
## Synthesizing constraints for a given deviation bound"""

# ╔═╡ 07135afb-b684-436b-824b-31a2a0d42ee0
md"""
Adjust the number of steps per tree $n$, and the total number of time steps $t$ for synthesizing the constraints.

|                     |                                                              |
|--------------------:|:-------------------------------------------------------------|
|               Model | $(@bind model_syn Select(model_names, default="ES"))         |
|                 $n$ | $(@bind n_syn Slider(1:20, default=10, show_value=true))     |
|                 $t$ | $(@bind time_syn Slider(5:5:500, default=100, show_value=true)) |
|     Max window size | $(@bind max_window_size_syn Slider(1:20, default=5, show_value=true))     |
|    Deviation bound  | $(@bind deviation_bound_syn Slider(5:50, default=10, show_value=true)) |
"""

# ╔═╡ 31f1cf8e-ee94-4b99-bd0a-f9ffc12924c6
function synthesize(deviation_bound, bounds, model, strat, n, max_window_size, time; dims=axes(model[1].A, 1))
	# get nominal behavior here
	
	deviations = fill(NaN, (max_window_size-1, max_window_size-1))
	# Windows size starting from 2

	min_hits = 1
	for window_size in 2:max_window_size
		# @info "Working on window_size: $(window_size)..."
		while min_hits < window_size
			automaton = strat(model[1], model[2], window_size-min_hits, window_size)
			augbounds = Augment(bounds, automaton)

			iterations = div(time+n-1, n)
			v, i = BoundedTreeIter(automaton, augbounds, n, iterations, deviation_bound, dims=dims)
			if v < deviation_bound
				# @info "  argmax(d) = $(i)"
				deviations[window_size-1, min_hits] = v
				break
			end

			@info "  maximum(d) = $(v)" window_size min_hits
			min_hits += 1
			
			# all_bounds = BoundedTreeIter(automaton, augbounds, n, iterations)
			# d = deviation(automaton, augbounds, all_bounds, dims=[1,2])
			# i = argmax(d)
			# v = maximum(d)

			# @info "  min_hits: $(min_hits) done."
			# if v < deviation_bound
			# 	deviations[window_size-1, min_hits] = v
			# 	break
			# end

			# min_hits += 1
		end
		# @info "window_size: $(window_size) done."
	end
	deviations
end

# ╔═╡ c4206822-0763-43d7-9d3a-ac3a200259ad
# ╠═╡ disabled = true
#=╠═╡
results = let
	q = size(model_map[model_syn][1].A, 1)
	bounds = repeat([10 10], q)
	synthesize(deviation_bound_syn, bounds, model_map[model_syn], HoldAndSkipNext, n_syn, max_window_size_syn, time_syn)
end
  ╠═╡ =#

# ╔═╡ 4173033e-d382-4c2f-8610-aee3436c7c28
function valid_weakly_hard_constraints(deviations)
	max_window_size = size(deviations, 1)
	results = []
	for min_hits in 1:max_window_size-1
		for window_size in reverse(min_hits+1:max_window_size)
			if !isnan(deviations[window_size-1, max_misses])
				push!(results, (min_hits, window_size))
				break
			end
		end
	end
end

# ╔═╡ 01c19e96-455b-4b51-a2d4-8dd6607254de
function schedule_two(ctrl_1, ctrl_2)
end

# ╔═╡ 46dd8536-65a6-4b31-84b2-efe231e462ac
md"""
## Examples
The reachable sets computed by the bounded tree method are shown below.  Try adjusting the parameters used for the model: the number of time steps per tree $n$, and the total number of time steps $t$.

|                     |                                                              |
|--------------------:|:-------------------------------------------------------------|
|               Model | $(@bind model_tree Select(model_names, default="ES"))        |
|                 $n$ | $(@bind n_trees Slider(1:20, default=10, show_value=true))   |
|          Max misses | $(@bind max_miss_trees Slider(-1:16, default=2, show_value=true)) |
|         Window size | $(@bind window_size_trees Slider(1:16, default=4, show_value=true)) |
|                 $t$ | $(@bind time_trees Slider(5:5:500, default=100, show_value=true)) |
| Sample trajectories | $(@bind sample_traj_trees CheckBox(default=false))           |
|    Deviation bound  | $(@bind deviation_bound_tree Slider(5:50, default=10, show_value=true)) |

First we generate the bounding boxes for possible trajectories.
"""

# ╔═╡ de8607d7-6eea-43de-acbe-ee3cb746f055
begin
	bounds = [1 1; 1 1]
	t_trees = div(time_trees+n_trees-1, n_trees)
	all_bounds = let
		model = model_map[model_tree]
		automaton = HoldAndSkipNext(model[1], model[2], max_miss_trees, window_size_trees)
		# automaton = HoldAndSkipNext(sysd, K, max_miss_trees, window_size_trees)
		augbounds = Augment(bounds, automaton)
		BoundedTreeIterOld(automaton, augbounds, n_trees, t_trees)
	end
end;

# ╔═╡ 6c091c25-9a81-4037-a9b7-37f0fdb899f1
let
	plot(title="Reachable set for $(time_trees) time steps, $(n_trees) per tree", xlabel="x_1", ylabel="x_2", legend=:bottomright)
	for t in 1:time_trees+1
		corners = corners_from_bounds(all_bounds[t,:,:], cycle=true)
		plot!(corners[1,:], corners[2,:], label=L"x[%$(t-1)]")
	end
	
	if sample_traj_trees
		model = model_map[model_tree]
		hsn = HoldAndSkipNext(model[1], model[2], 0, 1)
	    x, _ = Evol(Augment([bounds[1,1], bounds[2,1]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,1], bounds[2,2]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,2], bounds[2,1]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,2], bounds[2,2]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	end
	
	θ = LinRange(0, 2π, 500)
	plot!(legend=false)
end

# ╔═╡ 0d65cec3-6c42-4235-81ed-5781b6d5ea31
md"""
Then we compute the nominal behavior and obtain the deviation.
"""

# ╔═╡ 4a14195b-b93e-487a-b993-95dd8cff05e7
# ╠═╡ disabled = true
#=╠═╡
let
	model = model_map[model_tree]
	automaton = HoldAndSkipNext(model[1], model[2], max_miss_trees, window_size_trees)
	d = deviation(automaton, Augment(bounds, automaton), all_bounds, dims=[1, 2])
	i = argmax(d)
	v = maximum(d)
	
	plot(xlabel="Time", ylabel="Deviation", guidefontsize=18, tickfontsize=14, legendfontsize=14)
	if v < 1e3
		plot!(0:size(d,1)-1, d, label="Bounded Runs", linestyle=:dashdot, linewidth=2)
	else
		plot!(2:size(d,1)-1, d[3:end], label="Deviation", yaxis=:log)
	end
	scatter!([i-1], [v], label="", color=:black, series_annotations=[text(" ($(i-1), $(@sprintf("%.2f", v))) ", (i < n_trees*t_trees/2) ? :left : :right, 14)])
	plot!(legend=:right)
end
  ╠═╡ =#

# ╔═╡ 6ca31d15-a1d4-40bc-b534-35c870204bf4
md"""
## Appendix
"""

# ╔═╡ 06eea183-ffb2-4909-9d48-e28fe228b67e
md"""
### Miscellaneous Code
"""

# ╔═╡ b16aac7e-4b40-4f28-b9e1-d13ad64050f6
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ControlSystems = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[compat]
CSV = "~0.10.4"
ControlSystems = "~0.12.11"
Distances = "~0.10.7"
Distributions = "~0.25.56"
LaTeXStrings = "~1.3.0"
OffsetArrays = "~1.10.8"
Plots = "~1.28.1"
PlutoUI = "~0.7.38"
QuadGK = "~2.4.2"
Tables = "~1.7.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "c933ce606f6535a7c7b98e1d86d5d1014f730596"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "5.0.7"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "28bbdbf0354959db89358d1d79d421ff31ef0b5e"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.3"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "baaac45b4462b3b0be16726f38b789bf330fcb7a"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.21"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "f576084239e6bdf801007c80e27e2cc2cd963fe0"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.ControlSystems]]
deps = ["Colors", "DSP", "DelayDiffEq", "DiffEqCallbacks", "ForwardDiff", "IterTools", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MatrixEquations", "MatrixPencils", "OrdinaryDiffEq", "Polynomials", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "5c03162576ebf2d3f294960800ecded7c7c8a403"
uuid = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
version = "0.12.11"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "5e5f8f363c8c9a2415ef9185c4e0ff6966c87d52"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.2"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "3e03979d16275ed5d9078d50327332c546e24e68"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.5"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "UnPack"]
git-tree-sha1 = "52f54bd7f7bc1ce794add0ccf08f8fa21acfaed9"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.35.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "3c55535145325e0e3fa7a397e3a50e5f220d1edc"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.83.2"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "c4b99e3a199e293e7290eea94ba89364d47ee557"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.22.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "28d605d9a0ac17118fe2c5e9ce0fbb76c3ceb120"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "70f5bfdfbdc6c9d2b7a143d70ae88f4cb7b193b1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.56"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterface", "GenericSchur", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "789a5aaec6856ad774a95cdb558806afeeb256f9"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.15.0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "b6bf57ec7a3f294c97ae46124705a9e6b906a209"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.15"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "56956d1e4c1221000b7781104c58c34019792951"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "34e6147e7686a101c245f12dba43b743c7afda96"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.27"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "af237c08bda486b74318c8070adb96efa6952530"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cd6efcf9dc746b06709df14e462f0a3fe0786b1e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.2+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "18be5268cf415b5e27f34980ed25a7d34261aa83"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.7"

[[deps.Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "303d70c961317c4c20fafaf5dbe0e6d610c38542"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "cae5e3dfd89b209e01bcd65b3a25e74462c67ee0"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.3.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "82f5afb342a5624dc4651981584a841f6088166b"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.8.0"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "49b0c1dd5c292870577b8f58c51072bd558febb9"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.4"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "fbd884a02f8bf98fd90c53c1c9d2b21f9f30f42a"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.8.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "b651f573812d6c36c22c944dd66ef3ab2283dfa1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.6"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearMaps]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics"]
git-tree-sha1 = "1693d6d0dfefd24ee97ffc5ea91f1cd2cf77ef6e"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.6.1"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "6eb8e10ed29b85673495c29bd77ee0dfa8929977"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.15.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "76c987446e8d555677f064aaac1145c4c17662f8"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.14"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "4acc35e95bf18de5e9562d27735bef0950f2ed74"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.108"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MatrixEquations]]
deps = ["LinearAlgebra", "LinearMaps"]
git-tree-sha1 = "979680f9421531f4b8ff129849a108df8d395c17"
uuid = "99c1a7ee-ab34-5fd5-8076-27c950a045f4"
version = "2.2.1"

[[deps.MatrixPencils]]
deps = ["LinearAlgebra", "Polynomials", "Random"]
git-tree-sha1 = "864ae9033dc44114b112ee88752263cdd6a20f68"
uuid = "48965c70-4690-11ea-1f13-43a2532b2fa8"
version = "1.7.4"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "aeebff6a2a23506e5029fd2248a26aca98e477b3"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.16"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "8031a288c9b418664a3dfbac36e464a3f61ace73"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.10.0"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "3114946c67ef9925204cc024a73c9e679cebe0d7"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.8"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "d05baca9ec540de3d8b12ef660c7353aae9f9477"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.28.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "8d95a735921204f5d551ac300b20d802a150433a"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.8"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "7e597df97e46ffb1c8adbaddfa56908a7a20194b"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.5"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0107e2f7f90cc7f756fee8a304987c574bbd7583"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.0.0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "6c138c8510111fa47b5d2ed8ada482d97e279bee"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "bfe14f127f3e7def02a6c2b1940b39d0dabaa3ef"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.26.3"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "a9a852c7ebb08e2a40e8c0ab9830a744fa283690"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.10"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "ac399b5b163b9140f9c310dfe9e9aaa225617ff6"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.32"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "331f8250fe6a575522635cda5615a4e28fc60efe"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.31.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "314a07e191ea4a5ea5a2f9d6b39f03833bde5e08"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.21.0"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "91181e5820a400d1171db4382aa36e7fd19bee27"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "cd56bf18ed715e8b09f06ef8c6b781e6cdc49911"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c82aaa13b44ea00134f8c9c89819477bd3986ecd"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.3.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "70d9007ff05440058c0301985b2275edc2b2ce25"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.3.3"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "8f705dd141733d79aa2932143af6c6e0b6cea8df"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.6"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "b8d08f55b02625770c09615d96927b3a8396925e"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.11"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "858e541ffc21873e45aeaf744e0d015966e0328e"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.30"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═ee503656-cd73-11ec-13cb-61d08f418ed6
# ╠═a402d61b-d832-42b1-83f8-de95a87cd771
# ╟─39568986-8cf7-4fe2-81b9-8139f26f14fa
# ╠═110ca59b-f6e1-4033-9570-8e2b135f4f93
# ╟─3dc10fd0-9a52-4090-a3d5-4ff048c5a28b
# ╠═edd5649e-efcc-49bf-a333-0e14e8be518e
# ╟─d1f0cdfb-f6b6-4bac-a3e1-864be6e9be5e
# ╠═e2a404c9-ed92-4367-bad2-f425606b88dc
# ╟─8025bd64-5956-455b-b232-e41cb541ca09
# ╟─8e910c95-76dc-45c1-b09f-b0b725e49b7d
# ╟─ec833777-4b93-48e7-a076-dbcf147a908a
# ╠═41162c87-f54f-4e8c-8cfd-ac7a3916ec3d
# ╠═0216d280-02e6-4c24-b6e6-e1ff7b716e08
# ╟─fb97808c-0911-4d19-a29d-9ab97e9d3c5f
# ╠═fd65777a-0e7e-41a7-a11a-f6545bf9cef4
# ╟─083caada-cb9b-422f-9f1d-e72ce127f66e
# ╠═a6babf6a-48b0-4ff6-9b43-7b9871f506c8
# ╠═ce6cd63a-6393-4a3b-92c0-6ed3fb5edd79
# ╠═713a1e90-052d-42a7-b61d-8001da3c2129
# ╠═83347d7d-b559-4ee1-9272-2a249a378113
# ╠═de96915a-decc-4193-a63a-9964707dc025
# ╠═dc7b7a0e-3b0e-4c32-bc8e-9b1fe319702c
# ╟─839c0a6c-f6c2-4f53-abae-63739b4ad75b
# ╠═9ff9233d-d0b9-4687-a0d5-4a396a8bad54
# ╟─ec282dbc-e4e3-4197-8587-130b34ee9550
# ╠═805e54ef-be8a-4956-80a0-700ec1b96a92
# ╠═01dacf35-a7d4-4e2b-8391-d40fe3b317f7
# ╠═974d1e72-c6e8-4d61-acc8-a5263307bd51
# ╠═a43c4d8e-1935-4225-8db0-139b7b4207a6
# ╠═96a0379f-ffb3-47b4-8719-d60ada971d42
# ╠═b285273a-04fd-494f-bb15-ab0463e07496
# ╠═f6ebc370-9059-4236-836e-d080bd012fb9
# ╟─731f357a-a4c5-44dd-85c2-105a950d72e6
# ╠═811d1e05-a6ca-43bd-9f62-395392004747
# ╠═e014f895-6ed1-4489-8288-68e184f2ba39
# ╠═8c99c664-3cea-40a0-83c2-d8d1f6c8e46f
# ╟─48bfd3c8-87ba-482d-b0e2-22019e1afef4
# ╠═c2edd2a4-0ad1-45dc-8d6e-4c922a927e3c
# ╠═515802ce-3186-43d9-9897-6286622782c2
# ╟─a32a9411-1e3b-47dc-81f7-89ccf851d683
# ╠═550655ed-411c-4fbc-b932-730bc32514fc
# ╠═3873d023-15f1-490c-b3b7-caf6b96f8301
# ╟─80174b73-34f9-4c62-bf09-29672b68a558
# ╠═07135afb-b684-436b-824b-31a2a0d42ee0
# ╠═31f1cf8e-ee94-4b99-bd0a-f9ffc12924c6
# ╠═c4206822-0763-43d7-9d3a-ac3a200259ad
# ╠═4173033e-d382-4c2f-8610-aee3436c7c28
# ╠═01c19e96-455b-4b51-a2d4-8dd6607254de
# ╟─46dd8536-65a6-4b31-84b2-efe231e462ac
# ╠═de8607d7-6eea-43de-acbe-ee3cb746f055
# ╠═6c091c25-9a81-4037-a9b7-37f0fdb899f1
# ╠═0d65cec3-6c42-4235-81ed-5781b6d5ea31
# ╠═4a14195b-b93e-487a-b993-95dd8cff05e7
# ╠═6ca31d15-a1d4-40bc-b534-35c870204bf4
# ╟─06eea183-ffb2-4909-9d48-e28fe228b67e
# ╠═fd140d06-742c-44a9-87ee-304533065222
# ╠═b16aac7e-4b40-4f28-b9e1-d13ad64050f6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
