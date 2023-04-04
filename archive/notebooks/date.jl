### A Pluto.jl notebook ###
# v0.19.18

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

	using Plots, LaTeXStrings
	using PlutoUI
	using DelimitedFiles
	using DataStructures
	
	TableOfContents()
end

# ╔═╡ d9f35da9-b0ef-4faf-b550-9cd1aad448ad
md"""
# DATE 2023
Experimental results of our DATE 2023 paper: "Safety-aware Implementation of Control Tasks via Period Boosting and Compressing"

## Setting Up
"""

# ╔═╡ cb2f6baf-cd83-4d67-9041-2fae6badf40f
md"""
## Schedule Synthesis
"""

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

# ╔═╡ 7e7f7c54-316d-46d2-9520-0fb4fd440f28
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

# ╔═╡ 9b5e9ffa-e5d0-4313-b602-2be77f733ddd
"""
		digits_b2r(x[, pad])
	
	A shortcut for `digits(x, base=2, pad=pad) |> reverse`
	"""
	function digits_b2r(x::Int, pad::Int=0)
		digits(x, base=2, pad=pad) |> reverse
	end

# ╔═╡ aec6790b-50a2-4120-b01e-ee51b0355888
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
	
	function compare_gain(path1, path2, path3, threshold_percent, percentmode=true)
		dev1, ind1, w = getdev(path1, sigdigits=3)
		dev2, ind2, w = getdev(path2, sigdigits=3)
		dev3, ind3, w = getdev(path3, sigdigits=3)

		if percentmode
			min_v = minimum(x -> isnan(x) ? Inf  : x, vcat(dev1, dev2, dev3))
			max_v = maximum(x -> isnan(x) ? -Inf : x, vcat(dev1, dev2, dev3))
			threshold = min_v + (min(max_v, 1000) - min_v) * threshold_percent
		else
			threshold = threshold_percent
		end
		
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
		to_show
		
		# @info to_show
		# threshold
	end
	function compare_sys(path, sys, threshold, p1, p2)
		path1 = "$path/$(p1)s_$(sys["p"])s/"
		path2 = "$path/$(p2)s_$(sys["p"])s/"
		path3 = "$path/$(p1)s_$(p1)s/"
		filename = "$(sys["name"])_$(sys["x0"])_n$(sys["n"])_t100.csv"
		compare_gain(path1*filename, path2*filename, path3*filename, threshold, false)
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

# ╔═╡ b62022c8-c7f0-4ebe-a17b-a8eb43e796ca
begin
	t = 100
	
	RC_threshold = 1.4
	F1_threshold = 1.2
	DC_threshold = 3.5
	CS_threshold = 9.4
	CC_threshold = 0.53
end

# ╔═╡ ea5ec578-0d5e-4c40-8082-8b806d2b5e99
md"""
| System | WCET | Period | Safety Margin |
| ------ | ---- | ------ | ------------- |
| RC     | 10ms | 23ms   | $RC_threshold |
| F1     | 13ms | 20ms   | $F1_threshold |
| DCM    | 12ms | 23ms   | $DC_threshold |
| CSS    | 10ms | 27ms   | $CS_threshold |
| CC     | 15ms | 28ms   | $CC_threshold |
Utilization $$U \approx 2.51 > 1$$
"""

# ╔═╡ 24733f39-edb0-4222-b7d4-0235f885aeb7
systems = [
    Dict([("name", "RC"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", "delay_lqr"), ("ctrl_args", ()), ("th", RC_threshold)]),
    Dict([("name", "F1"), ("x0", 1), ("n", 16), ("p", 0.020), ("ctrl", "delay_lqr"), ("ctrl_args", ()), ("th", F1_threshold)]),
    Dict([("name", "DCM"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", "delay_lqr"), ("ctrl_args", ()), ("th", DC_threshold)]),
    Dict([("name", "CSS"), ("x0", 100), ("n", 15), ("p", 0.027), ("ctrl", "delay_lqr"), ("ctrl_args", ()), ("th", CS_threshold)]),
    Dict([("name", "CC2"), ("x0", 1), ("n", 20), ("p", 0.028), ("ctrl", "pole_place"), ("ctrl_args", (0.9)), ("th", CC_threshold)])
]

# ╔═╡ 9a60968f-31c1-4339-8d98-8d0dd0a16961
schedule_15_no, schedule_15_ya = let
	p1 = 0.015  # 40ms
	p2 = 0.028  # 28ms
	
	@info [compare_sys("../experiments/data/finalize", systems[1], RC_threshold, p1, p2), compare_sys("../experiments/data/finalize", systems[2], F1_threshold, p1, p2), compare_sys("../experiments/data/finalize", systems[3], DC_threshold, p1, p2), compare_sys("../experiments/data/ccs_lqr", systems[4], CS_threshold, p1, p2), compare_sys("../experiments/data/cc2_0.9pole", systems[5], CC_threshold, p1, p2)]

	RC1, RC2, RC3 = load_sys("../experiments/data/finalize", systems[1], RC_threshold, p1, p2)
	F11, F12, F13 = load_sys("../experiments/data/finalize", systems[2], F1_threshold, p1, p2)
	DC1, DC2, DC3 = load_sys("../experiments/data/finalize", systems[3], DC_threshold, p1, p2)
	CS1, CS2, CS3 = load_sys("../experiments/data/ccs_lqr", systems[4], CS_threshold, p1, p2)
	CC1, CC2, CC3 = load_sys("../experiments/data/cc2_0.9pole", systems[5], CC_threshold, p1, p2)
	
	schedule_15_no = nothing
	schedule_15_ya = nothing
	for c1 in RC1, c2 in F11, c3 in DC1, c4 in CS1, c5 in CC1
		path = isschedulable(c1, c2, c3, c4, c5, cores=1)
		if length(path) > 0
			schedule_15_no = map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			@info "Found results 15ms w/o redesign" c1 c2 c3 c4 c5 schedule_15_no
			break
		end
	end
	# for c1 in RC2, c2 in F12, c3 in DC2, c4 in CS2, c5 in CC2
	# 	path = isschedulable(c1, c2, c3, c4, c5, cores=2)
	# 	if length(path) > 0
	# 		@info "Found results 28ms w/o redesign" c1 c2 c3 c4 c5 map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
	# 		break
	# 	end
	# end
	for c1 in RC3, c2 in F13, c3 in DC3, c4 in CS3, c5 in CC3
		path = isschedulable(c1, c2, c3, c4, c5, cores=1)
		if length(path) > 0
			schedule_15_ya = map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			@info "Found results 15ms w/ redesign" c1 c2 c3 c4 c5 schedule_15_ya
			break
		end
	end
	"15ms and 28ms"
	schedule_15_no, schedule_15_ya
end

# ╔═╡ ccec61e4-5c71-49a4-bcd1-59655ba15311
schedule_28_no, schedule_28_ya = let
	p1 = 0.028  # 40ms
	p2 = 0.028  # 28ms
	
	@info [compare_sys("../experiments/data/finalize", systems[1], RC_threshold, p1, p2), compare_sys("../experiments/data/finalize", systems[2], F1_threshold, p1, p2), compare_sys("../experiments/data/finalize", systems[3], DC_threshold, p1, p2), compare_sys("../experiments/data/ccs_lqr", systems[4], CS_threshold, p1, p2), compare_sys("../experiments/data/cc2_0.9pole", systems[5], CC_threshold, p1, p2)]

	RC1, RC2, RC3 = load_sys("../experiments/data/finalize", systems[1], RC_threshold, p1, p2)
	F11, F12, F13 = load_sys("../experiments/data/finalize", systems[2], F1_threshold, p1, p2)
	DC1, DC2, DC3 = load_sys("../experiments/data/finalize", systems[3], DC_threshold, p1, p2)
	CS1, CS2, CS3 = load_sys("../experiments/data/ccs_lqr", systems[4], CS_threshold, p1, p2)
	CC1, CC2, CC3 = load_sys("../experiments/data/cc2_0.9pole", systems[5], CC_threshold, p1, p2)

	
	schedule_28_no = nothing
	schedule_28_ya = nothing
	
	for c1 in RC1, c2 in F11, c3 in DC1, c4 in CS1, c5 in CC1
		path = isschedulable(c1, c2, c3, c4, c5, cores=2)
		if length(path) > 0
			schedule_28_no = map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			@info "Found results 28ms w/o redesign" c1 c2 c3 c4 c5 schedule_28_no
			break
		end
	end
	for c1 in RC3, c2 in F13, c3 in DC3, c4 in CS3, c5 in CC3
		path = isschedulable(c1, c2, c3, c4, c5, cores=2)
		if length(path) > 0
			schedule_28_ya = map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			@info "Found results 28ms w/ redesign" c1 c2 c3 c4 c5 schedule_28_ya
			break
		end
	end
	"28ms"
	schedule_28_no, schedule_28_ya
end

# ╔═╡ d18c1d63-c647-41e3-a80f-7a0443365800
schedule_40_no, schedule_40_ya = let
	p1 = 0.040  # 40ms
	p2 = 0.028  # 28ms
	
	@info [compare_sys("../experiments/data/finalize", systems[1], RC_threshold, p1, p2), compare_sys("../experiments/data/finalize", systems[2], F1_threshold, p1, p2), compare_sys("../experiments/data/finalize", systems[3], DC_threshold, p1, p2), compare_sys("../experiments/data/ccs_lqr", systems[4], CS_threshold, p1, p2), compare_sys("../experiments/data/cc2_0.9pole", systems[5], CC_threshold, p1, p2)]

	RC1, RC2, RC3 = load_sys("../experiments/data/finalize", systems[1], RC_threshold, p1, p2)
	F11, F12, F13 = load_sys("../experiments/data/finalize", systems[2], F1_threshold, p1, p2)
	DC1, DC2, DC3 = load_sys("../experiments/data/finalize", systems[3], DC_threshold, p1, p2)
	CS1, CS2, CS3 = load_sys("../experiments/data/ccs_lqr", systems[4], CS_threshold, p1, p2)
	CC1, CC2, CC3 = load_sys("../experiments/data/cc2_0.9pole", systems[5], CC_threshold, p1, p2)
	
	schedule_40_no = nothing
	schedule_40_ya = nothing
	for c1 in RC1, c2 in F11, c3 in DC1, c4 in CS1, c5 in CC1
		path = isschedulable(c1, c2, c3, c4, c5, cores=3)
		if length(path) > 0
			schedule_40_no = map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			@info "Found results 40ms w/o redesign" c1 c2 c3 c4 c5 schedule_40_no
			break
		end
	end
	# for c1 in RC2, c2 in F12, c3 in DC2, c4 in CS2, c5 in CC2
	# 	path = isschedulable(c1, c2, c3, c4, c5, cores=2)
	# 	if length(path) > 0
	# 		schedule_28_no = map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
	# 		@info "Found results 28ms w/o redesign" c1 c2 c3 c4 c5 schedule_28_no
	# 		break
	# 	end
	# end
	for c1 in RC3, c2 in F13, c3 in DC3, c4 in CS3, c5 in CC3
		path = isschedulable(c1, c2, c3, c4, c5, cores=3)
		if length(path) > 0
			schedule_40_ya = map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
			@info "Found results 40ms w/ redesign" c1 c2 c3 c4 c5 schedule_40_ya
			break
		end
	end
	"40ms and 28ms"
	schedule_40_no, schedule_40_ya
end

# ╔═╡ bb040810-3ea8-4cf4-a5c4-252c156b16fc
md"""
## Diagrams for Dynamics
"""

# ╔═╡ efa8632b-c01a-45f3-a212-6d7f35c02045
"""
	ingredients(path::String)

Pluto workaround to import code: Instead of using `include()`, use this.
"""
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ bde5c5e5-e55c-4416-b453-cc8f7dd6cc06
function get_schedule(fullschedule, id, t)
	dupl = findfirst(x -> x==fullschedule[end], fullschedule)
	leng = length(fullschedule)
	ctrl_schedule = map(x -> x[id][end], fullschedule)
	# @info leng dupl

	res = zeros(Int64, t)
	for i = 1:t
		if i <= leng
			res[i] = ctrl_schedule[i]
		else
			res[i] = ctrl_schedule[((i-leng) % (leng-dupl)) + dupl]
		end
	end
	res
end

# ╔═╡ e50ebca1-84f3-4775-9ab6-a8cf48f1ac8c
exp = ingredients("../experiments/experiment.jl")

# ╔═╡ 44344f1f-9577-45d8-9f73-ab01e2ecff47
nominal = exp.nominal_traj(exp.sys_f1, exp.delay_lqr(exp.sys_f1, 0.02), 20, 1, 100)

# ╔═╡ 239a1e59-6793-446d-92c6-e1c214191178
md"""
|                   |                                                              |
|------------------:|:-------------------------------------------------------------|
|    Camera azimuth | $(@bind azimuth_2 Slider(-180:1:180, default=35, show_value=true)) |
|  Camera elevation | $(@bind elevation_2 Slider(-90:1:90, default=25, show_value=true)) |
"""

# ╔═╡ 5383ac12-0202-4fa7-90f7-45d3a9603891
function plotsts(systems, stsid, show, legend)
	sts = systems[stsid]
	
	hit_28_no = get_schedule(schedule_28_no, stsid, t-1) .+ 1
	hit_28_ya = get_schedule(schedule_28_ya, stsid, t-1) .+ 1
	# hit_40_no = get_schedule(schedule_40_no, stsid, t-1) .+ 1
	hit_40_ya = get_schedule(schedule_40_ya, stsid, t-1) .+ 1
	
	# traj_15_no = exp.create_traj(sts["name"], sts["x0"], t, (0.015, sts["p"]), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"], hitpattern=hit_28_no)
	# traj_15_ya = exp.create_traj(sts["name"], sts["x0"], t, (0.015, 0.015), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"], hitpattern=hit_28_ya)
	
	nominal = exp.create_traj(sts["name"], sts["x0"], t, (sts["p"], sts["p"]), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"])

	sys = exp.sys_map[sts["name"]]
	p_ms = Int64(sts["p"] * 1000)
	# nominal = exp.nominal_traj(sys, exp.delay_lqr(sys, sts["p"]), p_ms, sts["x0"], 100)
	@info nominal
	
	traj_28_no = exp.create_traj(sts["name"], sts["x0"], t, (0.028, sts["p"]), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"], hitpattern=hit_28_no)
	traj_28_ya = exp.create_traj(sts["name"], sts["x0"], t, (0.028, 0.028), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"], hitpattern=hit_28_ya)
	
	traj_40_no = exp.create_traj(sts["name"], sts["x0"], t, (0.040, sts["p"]), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"], hitpattern=2*ones(Int64, t-1))
	# traj_40_no = exp.create_traj(sts["name"], sts["x0"], t, (0.040, sts["p"]), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"], hitpattern=hit_40_ya)
	traj_40_ya = exp.create_traj(sts["name"], sts["x0"], t, (0.040, 0.040), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"], hitpattern=2*ones(Int64, t-1))

	traj_50_no = exp.create_traj(sts["name"], sts["x0"], t, (0.050, sts["p"]), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"])
	traj_50_ya = exp.create_traj(sts["name"], sts["x0"], t, (0.050, 0.050), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"])
	traj_60_no = exp.create_traj(sts["name"], sts["x0"], t, (0.060, sts["p"]), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"])
	traj_60_ya = exp.create_traj(sts["name"], sts["x0"], t, (0.060, 0.060), ctrl=sts["ctrl"], ctrl_args=sts["ctrl_args"])
	
	# plot(xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", legend=:topleft, format=:png, thickness_scaling=2, size=(1360,906))
	plt = plot(xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", legend=legend, camera=(azimuth_2,elevation_2), aspect_ratio=.7)
	
	# Plot safety pipe
	θ = LinRange(0, 2π, 40)
	circ_x = cos.(θ) * sts["th"]
	circ_y = sin.(θ) * sts["th"]
	for i in 1:t
		plot!(circ_x .+ nominal[i,1], circ_y .+ nominal[i,2], repeat([i], size(θ,1)), label=(i == 1) ? "Safety Margin" : "", seriestype=[:shape,], c=:lightblue, linecolor=:lightblue)
	end

	# # Plot random trajectories
	# label_good = false
	# # Plot the good trajectories first
	# for i = axes(trj_2, 1)
	# 	if dev_2[i] < pipe_radius
	# 		plot!(trj_2[i,:,1], trj_2[i,:,2], 1:H_2+1, label=(!label_good) ? "Random" : "", color=:green, opacity=10/n_samples_2)
	# 		label_good = true
	# 	end
	# end
	
	if show[1]
		plot!(traj_28_no[:,1], traj_28_no[:,2], (1:t).*0.028/sts["p"], label="28ms no redesign", color=:green, linewidth=2, opacity=0.7)
	end
	if show[2]
		plot!(traj_28_ya[:,1], traj_28_ya[:,2], (1:t).*0.028/sts["p"], label="28ms redesign", color=:green, linewidth=2, opacity=0.7)
	end
	if show[3]
		plot!(traj_40_no[:,1], traj_40_no[:,2], (1:t).*0.040/sts["p"], label="40ms no redesign", color=:red, linewidth=1.5, opacity=0.7)
	end
	if show[4]
		plot!(traj_40_no[1:20,1], traj_40_no[1:20,2], (1:20).*0.040/sts["p"], label="40ms no redesign", color=:red, linewidth=1.5, opacity=0.7)
	end
	if show[5]
		plot!(traj_40_ya[:,1], traj_40_ya[:,2], (1:t).*0.040/sts["p"], label="40ms redesign", color=:magenta, linewidth=2, opacity=0.7)
	end
	# plot!(traj_50_no[1:10,1], traj_50_no[1:10,2], 1:10, label="50 ms w/o redesign", color=:yellow, linewidth=1.5, opacity=0.7)
	# plot!(traj_50_ya[:,1], traj_50_ya[:,2], 1:t, label="50 ms w/ redesign", color=:pink, linewidth=2, opacity=0.7)
	
	# # Plot the bad trajectories on top of the good ones
	# label_bad = false
	nom_crosses = zeros(Float64, t)
	# for i = axes(trj_2, 1)
	# 	if dev_2[i] >= pipe_radius
	# 		t_viol = argmax(dist_2[i,:], dims=1)
	# 		#label_bad || plot!(cos.(θ) * pipe_radius .+ nominal_2[t_viol,1], sin.(θ) * pipe_radius .+ nominal_2[t_viol,2], repeat([t_viol,], size(θ,1)), label="Radius violated", style=:dot, linecolor=:gray)
	# 		crosses = [(d >= pipe_radius) ? 1 : 0 for d in dist_2[i,:]]
	# 		plot!(trj_2[i,:,1], trj_2[i,:,2], 1:H_2+1, label=(!label_bad) ? "Violation" : "", color=:red, markershape=:x, markeralpha=crosses)
	# 		nom_crosses .+= crosses
	# 		label_bad = true
	# 	end
	# end

	# # Finally, plot the nominal trajectory on top of everything else
	plot!(nominal[:,1], nominal[:,2], 1:t, label="Nominal Behavior", color=:black, linewidth=1, marker=:x, markeralpha=nom_crosses)

	plt
end

# ╔═╡ 802f13c0-a010-4386-8abd-70f6af2928cc
let
	# p1 = plotsts(systems, 2, [false, false, false, false, false], :topright)
	p2 = plotsts(systems, 5, [true, true, true, false, false], :topright)
	# plot(p1, p2, layout=2)
end

# ╔═╡ 36e98a22-9ab5-4662-bf72-f449824de4b8
p2 = plotsts(systems, 4, [false, false, false, false, false], :topright)

# ╔═╡ 68a80b94-ff3a-4c86-b997-7d295f58c889
md"""
## Old Experiments
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
# ╟─d9f35da9-b0ef-4faf-b550-9cd1aad448ad
# ╠═3b81315c-2a0b-11ed-0388-471cd2cb6d84
# ╟─cb2f6baf-cd83-4d67-9041-2fae6badf40f
# ╟─e1e49746-aa18-4064-83c8-cd2e773d38ee
# ╟─7e7f7c54-316d-46d2-9520-0fb4fd440f28
# ╟─9b5e9ffa-e5d0-4313-b602-2be77f733ddd
# ╟─aec6790b-50a2-4120-b01e-ee51b0355888
# ╟─58899f1b-03c9-4d0e-b49b-083b1d63b2fa
# ╟─84aca27f-c6d3-4824-892d-444fdc36d612
# ╟─f9c740f3-b134-455e-9472-eb091bda35e1
# ╟─ea5ec578-0d5e-4c40-8082-8b806d2b5e99
# ╠═24733f39-edb0-4222-b7d4-0235f885aeb7
# ╠═b62022c8-c7f0-4ebe-a17b-a8eb43e796ca
# ╟─9a60968f-31c1-4339-8d98-8d0dd0a16961
# ╟─ccec61e4-5c71-49a4-bcd1-59655ba15311
# ╠═d18c1d63-c647-41e3-a80f-7a0443365800
# ╟─bb040810-3ea8-4cf4-a5c4-252c156b16fc
# ╟─efa8632b-c01a-45f3-a212-6d7f35c02045
# ╠═bde5c5e5-e55c-4416-b453-cc8f7dd6cc06
# ╠═e50ebca1-84f3-4775-9ab6-a8cf48f1ac8c
# ╠═44344f1f-9577-45d8-9f73-ab01e2ecff47
# ╟─239a1e59-6793-446d-92c6-e1c214191178
# ╠═802f13c0-a010-4386-8abd-70f6af2928cc
# ╠═36e98a22-9ab5-4662-bf72-f449824de4b8
# ╠═5383ac12-0202-4fa7-90f7-45d3a9603891
# ╟─68a80b94-ff3a-4c86-b997-7d295f58c889
# ╟─470a7307-82a4-4a30-a308-6d8ade16a2f4
# ╠═57c2fd32-48b2-4a53-9bb3-96f28ef6c70b
# ╟─f17e5f30-8dbb-4079-a22b-d79b05aedbed
# ╠═421842ca-442c-499f-8ff5-1ac46d53e67e
# ╟─e3f33864-27e5-4bff-be7f-942d3a8cdaf8
# ╟─f37293d5-3ba0-4abf-898a-8d500d11302e
# ╟─4d862328-ca2f-4656-a296-e08237fd3150
# ╟─b3dd90f5-c87d-4d67-8094-3eaf836033f0
# ╟─076dc76e-b93c-4642-a189-5bf233ba4eb0
# ╟─adb5ed1f-d1ab-4eda-a362-5fd26e98b6c7
# ╟─71875c91-2167-4fa6-a26d-91f779269796
# ╟─29e7a702-5ef4-4f1f-bfcb-5a0b69621749
