### A Pluto.jl notebook ###
# v0.19.12

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
			# @info "Utilization is ok but no schedule found" constraints current_states
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

# ╔═╡ f9c740f3-b134-455e-9472-eb091bda35e1
begin
	function getdev(path; sigdigits=2)
		d = readdlm(path, ',', Float64)
		d = round.(d; sigdigits=sigdigits)
		w = size(d, 2)
		dev = d[1:w,:]
		ind = d[w+1:2w,:]
	
		dev, ind, w
	end
	
	"""
		load_result(path, threashold)
	
	Load the result (in `.csv` format) from the file system, and present the result according to the threashold value.
	"""
	function load_result(path, threshold_ratio, percentmode=true)
		deviations, indices, w = getdev(path, sigdigits=3)
		# @info deviations

		if percentmode
			min_v = minimum(x -> isnan(x) ? Inf  : x, deviations)
			max_v = maximum(x -> isnan(x) ? -Inf : x, deviations)
			threshold = min_v + (min(max_v, 1000) - min_v) * threshold_ratio
		else
			threshold = threshold_ratio
		end
	
		to_show = fill(NaN, size(deviations))
		constraints = Vector{WeaklyHardConstraint}()
		
		min_hits = 1
		for window_size in 2:w+1
			while min_hits < window_size
				if deviations[window_size-1, min_hits] <= threshold && 
				   indices[window_size-1, min_hits] <= 80
					to_show[window_size-1, min_hits] = deviations[window_size-1, min_hits]
					push!(constraints, WeaklyHardConstraint(min_hits, window_size))
					break
				end
				min_hits += 1
			end
		end
		
		# @info "Full results" to_show deviations indices
		# @info "Stats" deviations threshold
		constraints
	end
	
	function compare_gain(path1, path2, path3, threshold_percent)
		dev1, ind1, w = getdev(path1, sigdigits=3)
		dev2, ind2, w = getdev(path2, sigdigits=3)
		dev3, ind3, w = getdev(path3, sigdigits=3)
	
		min_v = minimum(x -> isnan(x) ? Inf  : x, vcat(dev1, dev2, dev3))
		max_v = maximum(x -> isnan(x) ? -Inf : x, vcat(dev1, dev2, dev3))
		threshold = min_v + (min(max_v, 1000) - min_v) * threshold_percent
		
		to_show = fill("   ", (w, w))
	
		min_hits = 1
		for window_size in 2:w+1
			for min_hits in 1:window_size-1
				d1, i1 = dev1[window_size-1, min_hits], ind1[window_size-1, min_hits]
				d2, i2 = dev2[window_size-1, min_hits], ind2[window_size-1, min_hits]
				d3, i3 = dev3[window_size-1, min_hits], ind3[window_size-1, min_hits]
				sy = ""
				sy *= if d1 <= threshold && i1 <= 80 "1" else " " end
				sy *= if d2 <= threshold && i2 <= 80 "2" else " " end
				sy *= if d3 <= threshold && i3 <= 80 "3" else " " end
				to_show[window_size-1, min_hits] = sy
				# min_hits += 1
			end
		end
		
		# @info deviations
		# @info deviations indices
		# display("text/plain", [dev1; ind1])
		# @info threshold
		@info to_show
		threshold
	end
	function compare_sys(path, sys, threshold, p1, p2)
		path1 = "$path/$(p1)s_$(sys["p"])s/"
		path2 = "$path/$(p2)s_$(sys["p"])s/"
		path3 = "$path/$(p1)s_$(p1)s/"
		filename = "$(sys["name"])_$(sys["x0"])_n$(sys["n"])_t100.csv"
		compare_gain(path1*filename, path2*filename, path3*filename, threshold)
	end
	function load_sys(path, sys, threshold, p1, p2)
		path1 = "$path/$(p1)s_$(sys["p"])s/"
		path2 = "$path/$(p2)s_$(sys["p"])s/"
		path3 = "$path/$(p1)s_$(p1)s/"
		filename = "$(sys["name"])_$(sys["x0"])_n$(sys["n"])_t100.csv"
		load_result(path1*filename, threshold, false), 
		load_result(path2*filename, threshold, false),
		load_result(path3*filename, threshold, false)
	end

	[getdev, load_result, compare_gain, compare_sys, load_sys]
end

# ╔═╡ 24733f39-edb0-4222-b7d4-0235f885aeb7
systems = [
    Dict([("name", "RC"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    Dict([("name", "F1"), ("x0", 1), ("n", 16), ("p", 0.020), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    Dict([("name", "DCM"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    Dict([("name", "CSS"), ("x0", 100), ("n", 15), ("p", 0.027), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    Dict([("name", "CC2"), ("x0", 1), ("n", 20), ("p", 0.028), ("ctrl", "pole_place"), ("ctrl_args", (0.85))])
]

# ╔═╡ 1a8a6ddb-c4c6-4ae4-98c8-870d4ce8448c
begin
	p1 = 0.040  # 40ms
	p2 = 0.028  # 28ms

	# case 1: period = 40ms (3 tasks), w/o redesign
	# case 2: period = 28ms, w/o redesign
	# case 3: period = 40ms, w/  redesign (mention F1 is the main culprit if have space)
	
	# case 4: period = 28ms, w/  redesign
	# case 5: period = 15ms, w/o redesign
	# case 6: period = 15ms, w/  redesign
end

# ╔═╡ ffd30614-bcd1-4731-95c8-53a8afeb42da
@bind RC Slider(0:0.01:1, default=0.32, show_value=true)

# ╔═╡ 1f6cb76b-8d2b-49a1-a8b6-3d16e617a2c0
RC_threshold = compare_sys("../experiments/data/finalize", systems[1], RC, p1, p2)

# ╔═╡ 1ae579f7-4305-434d-815a-1a2fe693d7ce
@bind F1 Slider(0:0.0001:0.001, default=0.02, show_value=true)

# ╔═╡ 0aa39f38-1c75-45aa-bb24-1aaa47128535
F1_threshold = compare_sys("../experiments/data/finalize", systems[2], F1, p1, p2)

# ╔═╡ 2df35efa-7df9-4ec8-b2f7-9c0b3828ad1e
@bind DC Slider(0:0.01:1, default=0.5, show_value=true)

# ╔═╡ 2f15b86d-193e-45fa-a99b-b10f118f69fe
DC_threshold = compare_sys("../experiments/data/finalize", systems[3], DC, p1, p2)

# ╔═╡ 909d3c26-1030-42f8-b46b-d3cb4b1e7a70
@bind CS Slider(0:0.01:10, default=0.76, show_value=true)

# ╔═╡ 641a246c-0fdb-4c5d-8d6a-5c17617b6904
CS_threshold = compare_sys("../experiments/data/ccs_lqr", systems[4], CS, p1, p2)

# ╔═╡ 463edf21-d243-4337-9475-528cd3a1bb1a
@bind CC Slider(0:0.0001:0.001, default=0.0004, show_value=true)

# ╔═╡ 0e44d745-a73b-4ff5-b38d-7f94d154d508
CC_threshold = compare_sys("../experiments/data/cc2_0.9pole", systems[5], CC, p1, p2)

# ╔═╡ 5f18424a-5085-4e43-a3ba-6f7e51a3e2dc
let
	RC1, RC2, RC3 = load_sys("../experiments/data/finalize", systems[1], RC_threshold, p1, p2)
	F11, F12, F13 = load_sys("../experiments/data/finalize", systems[2], F1_threshold, p1, p2)
	DC1, DC2, DC3 = load_sys("../experiments/data/finalize", systems[3], DC_threshold, p1, p2)
	CS1, CS2, CS3 = load_sys("../experiments/data/ccs_lqr", systems[4], CS_threshold, p1, p2)
	CC1, CC2, CC3 = load_sys("../experiments/data/cc2_0.9pole", systems[5], CC_threshold, p1, p2)

	# @info RC1, F11, DC1, CS1, CC1
	# @info RC2, F12, DC2, CS2, CC2
	# @info RC3, F13, DC3, CS3, CC3
	
	for c1 in RC1, c2 in F11, c3 in DC1, c4 in CS1, c5 in CC1
		path = isschedulable(c1, c2, c3, c4, c5, cores=3)
		if length(path) > 0
			@info "Found results case 1" c1 c2 c3 c4 c5 map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			break
		end
	end
	for c1 in RC2, c2 in F12, c3 in DC2, c4 in CS2, c5 in CC2
		path = isschedulable(c1, c2, c3, c4, c5, cores=2)
		if length(path) > 0
			@info "Found results case 2" c1 c2 c3 c4 c5 map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			break
		end
	end
	for c1 in RC3, c2 in F13, c3 in DC3, c4 in CS3, c5 in CC3
		path = isschedulable(c1, c2, c3, c4, c5, cores=3)
		if length(path) > 0
			@info "Found results case 3" c1 c2 c3 c4 c5 map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			break
		end
	end
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
prefix = "../experiments/data/archive/common_period/"

# ╔═╡ 57c2fd32-48b2-4a53-9bb3-96f28ef6c70b
@bind threshold_rc Slider(0:0.01:1, default=0.58, show_value=true)

# ╔═╡ f17e5f30-8dbb-4079-a22b-d79b05aedbed
rc = load_result("../experiments/data/archive/common_period_names/0.028s_0.023s/RC_100_n5_t100.csv", threshold_rc)

# ╔═╡ 421842ca-442c-499f-8ff5-1ac46d53e67e
@bind threshold_f1 Slider(0:0.01:4, default=3.41, show_value=true)

# ╔═╡ e3f33864-27e5-4bff-be7f-942d3a8cdaf8
f1 = load_result("../experiments/data/archive/common_period_names/0.028s_0.02s/F1_1_n10_t100.csv", threshold_f1, false)

# ╔═╡ f37293d5-3ba0-4abf-898a-8d500d11302e
@bind threshold_dcm Slider(0:0.01:1, default=0.66, show_value=true)

# ╔═╡ 4d862328-ca2f-4656-a296-e08237fd3150
dcm = load_result("../experiments/data/archive/common_period_names/0.028s_0.023s/DCM_100_n5_t100.csv", threshold_dcm)

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
# ╠═e1e49746-aa18-4064-83c8-cd2e773d38ee
# ╠═9e4e98f9-ea74-40a5-8fad-ff90002a6a27
# ╠═58899f1b-03c9-4d0e-b49b-083b1d63b2fa
# ╠═84aca27f-c6d3-4824-892d-444fdc36d612
# ╠═f9c740f3-b134-455e-9472-eb091bda35e1
# ╠═24733f39-edb0-4222-b7d4-0235f885aeb7
# ╠═1a8a6ddb-c4c6-4ae4-98c8-870d4ce8448c
# ╠═ffd30614-bcd1-4731-95c8-53a8afeb42da
# ╠═1f6cb76b-8d2b-49a1-a8b6-3d16e617a2c0
# ╠═1ae579f7-4305-434d-815a-1a2fe693d7ce
# ╠═0aa39f38-1c75-45aa-bb24-1aaa47128535
# ╠═2df35efa-7df9-4ec8-b2f7-9c0b3828ad1e
# ╠═2f15b86d-193e-45fa-a99b-b10f118f69fe
# ╠═909d3c26-1030-42f8-b46b-d3cb4b1e7a70
# ╠═641a246c-0fdb-4c5d-8d6a-5c17617b6904
# ╠═463edf21-d243-4337-9475-528cd3a1bb1a
# ╠═0e44d745-a73b-4ff5-b38d-7f94d154d508
# ╠═5f18424a-5085-4e43-a3ba-6f7e51a3e2dc
# ╟─ea5ec578-0d5e-4c40-8082-8b806d2b5e99
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
