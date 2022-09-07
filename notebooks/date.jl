### A Pluto.jl notebook ###
# v0.19.11

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

# ╔═╡ 3b81315c-2a0b-11ed-0388-471cd2cb6d84
begin
	import Pkg
	Pkg.activate("..")
	
	using PlutoUI
	using DelimitedFiles
	using DataStructures
end

# ╔═╡ d2fe79ad-3b5f-461c-be6f-f3e7aaafe8fa
# TableOfContents()

# ╔═╡ e1e49746-aa18-4064-83c8-cd2e773d38ee
begin
	struct WeaklyHardConstraint
		minhits::Int 		# Minimum # of hits in a window
		windowsize::Int 	# The size of the window
	end
	struct Automaton
		L::Int64			# # of locations. Legal locations are in the range 0:L-1.
		Σ::Int64			# # of actions. Legal actions are in the range 0:Σ-1.
		T::Function			# Transition function. T(l,σ) is a location in 0:L-1.
		l_0::Int64			# Initial location in L.
		Q_f::Function		# Function that returns whether a location is final
	end
	struct SynthesizedAutomaton
		N::Int64			# # of comprising automata
		B::Vector{Int64}	# List of bits for all comprising controllers
		L::Int64			# Locations. Legal locations are in the range 0:L-1.
		Σ::Vector{Int64}	# List of actions.
		T::Function			# Transition function.  T(l,σ) is a location in 0:L-1.
		l_0::Int64			# Initial location in L.
		Q_f::Function		# Function that returns whether a location is final.
	end

	[WeaklyHardConstraint, Automaton, SynthesizedAutomaton]
end

# ╔═╡ 9e4e98f9-ea74-40a5-8fad-ff90002a6a27
begin
	"""
		undigit(d[, base=2])
	
	Convert a list of digits to a number. Default base=2.
	```julia
	julia> undigit([1, 0, 0])
	4
	```
	"""
	function undigit(d; base=2)
	    s = zero(eltype(d))
	    mult = one(eltype(d))
	    for val in reverse(d)
	        s += val * mult
	        mult *= base
	    end
	    return s
	end

	"""
		digits_b2r(x[, pad])
	
	A shortcut for `digits(x, base=2, pad=pad) |> reverse`
	"""
	function digits_b2r(x::Int, pad::Int=0)
		digits(x, base=2, pad=pad) |> reverse
	end
	
	"""
		state_separation(l, B[, indigits=false])
	
	state_separation takes a number `l` representing the overall state of
	a `SynthesizedAutomaton` and an array `B` representing the number of bits
	for each comprising `Automaton` of the `SynthesizedAutomaton` to compute a 
	list of individual states. For example
	
	```julia
	julia> state_separation(6, [1, 2])
	[1, 2]
	```
	The explanation of the example is as follows
	```
	6 -> [1, 1, 0]
	[1, 1, 0] / [1, 2] -> [[1], [1, 0]]
	[[1], [1, 0]] -> [1, 2]
	```
	"""
	function state_separation(l, B; indigits=false)
		@assert l < 2^sum(B) "l has more bits than the sum of B"
		
		bits = digits(l, base=2, pad=sum(B)) |> reverse
		index_points = cumsum(B)
		bits_separated = [bits[i+1:j] for (i, j) in zip([0; index_points[1:end-1]], index_points)]
		if indigits
			state_separated = bits_separated
		else
			state_separated = map(undigit, bits_separated)
		end
	end

	[undigit, digits_b2r, state_separation]
end

# ╔═╡ 58899f1b-03c9-4d0e-b49b-083b1d63b2fa
begin
	function AbstractControllerAutomaton(c::WeaklyHardConstraint)
		# Define the automaton's structure. State l=0 is the dead state (trapping)
		if c.minhits == 0					# No requirement. Always valid.
			L = 2
			T = (l, σ) -> 1
		elseif c.minhits == c.windowsize	# No misses allowed.
			L = 2
			T = (l, σ) -> l & σ
		else								# Dead state is l = 1
			L = 2^c.windowsize
			T = function (l, σ)
				l_new = (l << 1) & (L - 1) | σ
				l == 0 || count_ones(l_new) < c.minhits ? 0 : l_new
			end
		end
		
		# Put it all together. Starting position is L-1 since there is no miss from the beginning
		Automaton(L, 2, T, L-1, !iszero)
	end
	
	function SchedulerAutomaton(controllers::Automaton...; cores=1)
		# Converting tuple to array with collect()
		N = length(controllers)
		B = map(a -> ceil(Int64, log2(a.L)), controllers) |> collect
		L = 2^sum(B)
	
		all_actions = 0:2^length(controllers)-1
		Σ = filter(σ -> count_ones(σ) <= cores, all_actions)
	
		function T(l::Int, σ::Int)
			@assert l < L "Illegal location"
			@assert σ in Σ "Illegal action"
			states = state_separation(l, B)
			actions = digits_b2r(σ, N)
			# @info σ, actions
			new_locations = map((controller, l, σ) -> controller.T(l, σ), controllers, states, actions)
			# @info states, actions, new_locations
			new_location_bits = map((x, y) -> digits_b2r(x, y), new_locations, B)
			# @info new_location_bits
			result = Iterators.flatten(new_location_bits) |> collect |> undigit
			# @info digits_b2r(result)
			result
		end
	
		function Q_f(l::Int)
			states = state_separation(l, B);
			all(map((c, l) -> c.Q_f(l), controllers, states))
		end
	
		SynthesizedAutomaton(N, B, L, Σ, T, L-1, Q_f)
	end

	[AbstractControllerAutomaton, SchedulerAutomaton]
end

# ╔═╡ 84aca27f-c6d3-4824-892d-444fdc36d612
begin
	"""
		isschedulable(constraints...; cores=1, steps=100)

	Test if a set of tasks are schedulable. Each task has a list of constraints and and at each time slot, at most a few number of controllers can be run.
	"""
	function isschedulable(constraints...; cores=1, steps=100)
		utilization = sum(c -> c.minhits/c.windowsize, constraints)
		if utilization > cores
			# @info "Utilization exceeding number of cores" utilization
			return Vector{Int64}()
		end
	
		controllers = AbstractControllerAutomaton.(constraints)
		a = SchedulerAutomaton(controllers..., cores=cores)
	
		# cycle_detect = Set{Int64}(a.l_0)
		# path = Dict{Int64, Vector{Int64}}(a.l_0=>[a.l_0])
		current_states = Dict{Int64, Cons{Int64}}(a.l_0=>list(a.l_0))
	
		for step in 1:steps
			next_states = Dict{Int64, Cons{Int64}}()
			for (l, path) in pairs(current_states), σ in a.Σ
				l_new = a.T(l, σ)
				if !a.Q_f(l_new) || haskey(next_states, l_new)
					continue
				elseif l_new in path
					# @info "Cycle detected, returning early"
					return reverse(cons(l_new, path))
				end
				next_states[l_new] = cons(l_new, path)
			end
			if isempty(next_states)
				# some_non_readable_nonsense = [map(l -> state_separation(l, [c.windowsize for c in constraints], indigits=true), reverse(path)) |> collect for path in values(current_states)]
				# @info "Deat at step=$step" constraints some_non_readable_nonsense[1] some_non_readable_nonsense[2] some_non_readable_nonsense[3] some_non_readable_nonsense[4]
				@info "Utilization is ok but no schedule found" constraints current_states
				return Vector{Int64}()
			end
				
			current_states = next_states
		end
	
		for (l, path) in pairs(current_states)
			if a.Q_f(l)
				return reverse(path)
			end
		end
	
		@info "Utilization is ok but no schedule found" constraints current_states
		Vector{Int64}()
	end
	
	"""
		load_result(path, threashold)
	
	Load the result (in `.csv` format) from the file system, and present the result according to the threashold value.
	"""
	function load_result(path, threshold_ratio, percentmode=true)
		d = readdlm(path, ',', Float64)
		d = round.(d; sigdigits=3)
		w = size(d, 2)
		deviations = d[1:w,:]
		indices    = d[w+1:2w,:]
		# @info deviations

		if percentmode
			min_v = minimum(x -> isnan(x) ? Inf  : x, deviations)
			max_v = maximum(x -> isnan(x) ? -Inf : x, deviations)
			threshold = round(min_v + (max_v - min_v) * threshold_ratio, sigdigits=3)
		else
			threshold = threshold_ratio
		end
	
		to_show = fill(NaN, size(deviations))
		constraints = Vector{WeaklyHardConstraint}()
		
		min_hits = 1
		for window_size in 2:w+1
			while min_hits < window_size
				if deviations[window_size-1, min_hits] <= threshold
					to_show[window_size-1, min_hits] = deviations[window_size-1, min_hits]
					push!(constraints, WeaklyHardConstraint(min_hits, window_size))
					break
				end
				min_hits += 1
			end
		end
		
		# @info "Full results" to_show deviations indices
		@info "Stats" deviations threshold
		constraints
	end

	[isschedulable, load_result]
end

# ╔═╡ ea5ec578-0d5e-4c40-8082-8b806d2b5e99
md"""
| System | WCET | Period | Safety Margin |
| ------ | ---- | ------ | ------------- |
| RC     | 10ms | 23ms   | 1.59          |
| F1     | 13ms | 20ms   | 3.41          |
| DCM    | 12ms | 23ms   | 3.36          |
| CSS    | 10ms | 27ms   | 57.5          |
| CC     | 15ms | 28ms   | 1.05          |
Utilization $$U \approx 2.51 > 1$$
"""

# ╔═╡ 470a7307-82a4-4a30-a308-6d8ade16a2f4
prefix = "../experiments/data/common_period/"

# ╔═╡ 57c2fd32-48b2-4a53-9bb3-96f28ef6c70b
@bind threshold_rc Slider(0:0.01:1, default=0.58, show_value=true)

# ╔═╡ f17e5f30-8dbb-4079-a22b-d79b05aedbed
rc = load_result("../experiments/data/common_period_names/0.028s_0.023s/RC_100_n5_t100.csv", threshold_rc)

# ╔═╡ 421842ca-442c-499f-8ff5-1ac46d53e67e
@bind threshold_f1 Slider(0:0.01:4, default=3.41, show_value=true)

# ╔═╡ e3f33864-27e5-4bff-be7f-942d3a8cdaf8
f1 = load_result("../experiments/data/common_period_names/0.028s_0.02s/F1_1_n10_t100.csv", threshold_f1, false)

# ╔═╡ f37293d5-3ba0-4abf-898a-8d500d11302e
@bind threshold_dcm Slider(0:0.01:1, default=0.66, show_value=true)

# ╔═╡ 4d862328-ca2f-4656-a296-e08237fd3150
dcm = load_result("../experiments/data/common_period_names/0.028s_0.023s/DCM_100_n5_t100.csv", threshold_dcm)

# ╔═╡ b3dd90f5-c87d-4d67-8094-3eaf836033f0
@bind threshold_css Slider(0:0.01:1, default=0.34, show_value=true)

# ╔═╡ 076dc76e-b93c-4642-a189-5bf233ba4eb0
css = load_result(prefix * "0.028s_0.027s/CC2_HoldAndKill_n5_t100.csv", threshold_css)

# ╔═╡ adb5ed1f-d1ab-4eda-a362-5fd26e98b6c7
@bind threshold_cc2 Slider(0:0.01:1, default=0.0, show_value=true)

# ╔═╡ 71875c91-2167-4fa6-a26d-91f779269796
cc2 = load_result(prefix * "0.028s_0.028s/CC2_HoldAndKill_n5_t100.csv", threshold_cc2)

# ╔═╡ 29e7a702-5ef4-4f1f-bfcb-5a0b69621749
for c1 in rc, c2 in f1, c3 in dcm, c4 in css, c5 in cc2
	path = isschedulable(c1, c2, c3, c4, c5, cores=2)
	if length(path) > 0
		@info "Found results" c1 c2 c3 c4 c5 map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
		break
	end
end

# ╔═╡ Cell order:
# ╠═3b81315c-2a0b-11ed-0388-471cd2cb6d84
# ╠═d2fe79ad-3b5f-461c-be6f-f3e7aaafe8fa
# ╟─e1e49746-aa18-4064-83c8-cd2e773d38ee
# ╟─9e4e98f9-ea74-40a5-8fad-ff90002a6a27
# ╟─58899f1b-03c9-4d0e-b49b-083b1d63b2fa
# ╠═84aca27f-c6d3-4824-892d-444fdc36d612
# ╠═ea5ec578-0d5e-4c40-8082-8b806d2b5e99
# ╠═470a7307-82a4-4a30-a308-6d8ade16a2f4
# ╠═57c2fd32-48b2-4a53-9bb3-96f28ef6c70b
# ╠═f17e5f30-8dbb-4079-a22b-d79b05aedbed
# ╠═421842ca-442c-499f-8ff5-1ac46d53e67e
# ╠═e3f33864-27e5-4bff-be7f-942d3a8cdaf8
# ╠═f37293d5-3ba0-4abf-898a-8d500d11302e
# ╠═4d862328-ca2f-4656-a296-e08237fd3150
# ╠═b3dd90f5-c87d-4d67-8094-3eaf836033f0
# ╠═076dc76e-b93c-4642-a189-5bf233ba4eb0
# ╠═adb5ed1f-d1ab-4eda-a362-5fd26e98b6c7
# ╠═71875c91-2167-4fa6-a26d-91f779269796
# ╠═29e7a702-5ef4-4f1f-bfcb-5a0b69621749
