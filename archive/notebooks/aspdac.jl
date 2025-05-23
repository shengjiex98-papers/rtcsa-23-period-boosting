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

# ╔═╡ 3fbc35e1-e930-4a60-9eb8-283d733901d3
begin
	import Pkg
	Pkg.activate("..")
	
	using PlutoUI
	using DataStructures
end

# ╔═╡ b1c60c0a-e804-461b-a026-f173d963599c
begin
	using DelimitedFiles
end

# ╔═╡ 191d8766-d7d4-11ec-2c79-c5eee06d2d64
md"""
# Scheduler Automata
Setting up the notebook
"""

# ╔═╡ b6f14d07-77aa-4948-af8e-7d6062d6c38f
TableOfContents()

# ╔═╡ df7c19a5-69a4-4bd6-a157-572829150052
md"""
## Defining automata

We define two automata: first one corresponds to the regular language that describes a valid constraint for a single controller, while the second one describes a scheduler automaton that will be automatically generated given individual automata for different controllers.
"""

# ╔═╡ 2c3b7ca2-cda4-4274-a14c-8a57c16d15cd
struct WeaklyHardConstraint
	minhits::Int 		# Minimum # of hits in a window
	windowsize::Int 	# The size of the window
end

# ╔═╡ 3af3641c-687b-42bb-bfea-0f73425cb827
struct Automaton
	L::Int64			# # of locations. Legal locations are in the range 0:L-1.
	Σ::Int64			# # of actions. Legal actions are in the range 0:Σ-1.
	T::Function			# Transition function. T(l,σ) is a location in 0:L-1.
	l_0::Int64			# Initial location in L.
	Q_f::Function		# Function that returns whether a location is final
end

# ╔═╡ f1de5c80-c5bf-4d8d-bbb1-78ce7d675ff5
md"""
Next, we define the synthesized automaton. It takes as input a number of controllers and optionally the number of controllers that can be scheduled simultaneously (# of cores) and output a synthesized scheduler automaton.
"""

# ╔═╡ e7bf324a-c672-4161-84c9-4993efd8a0af
struct SynthesizedAutomaton
	N::Int64			# # of comprising automata
	B::Vector{Int64}	# List of bits for all comprising controllers
	L::Int64			# Locations. Legal locations are in the range 0:L-1.
	Σ::Vector{Int64}	# List of actions.
	T::Function			# Transition function.  T(l,σ) is a location in 0:L-1.
	l_0::Int64			# Initial location in L.
	Q_f::Function		# Function that returns whether a location is final.
end

# ╔═╡ b39b4969-03df-41d9-99b6-fb95d98a08ba
md"""
## Defining abstract controller automaton

Before we define the automaton for abstract controller, let's define a helper function that checks a hit/miss pattern against a constraint.
"""

# ╔═╡ ce56b597-02c5-4d3d-9712-762b493ae7a4
# ╠═╡ disabled = true
#=╠═╡
function check_constraint(x::Int64, constraint::WeaklyHardConstraint)
	mask = 2^constraint.windowsize - 1
	count_ones(x & mask) >= constraint.minhits
end
  ╠═╡ =#

# ╔═╡ d8f5a9cf-d535-40c4-b193-54d1b8a00e50
# ╠═╡ disabled = true
#=╠═╡
function check_constraint(x::Int64, constraints::Vector{WeaklyHardConstraint})
	all(check_constraint.(x, constraints))
end
  ╠═╡ =#

# ╔═╡ 5da180a4-7a0e-463d-bbb3-74e83a23eb47
md"""
Note that the below example returns true. This is because `count_ones_window` only considers the newest `windowsize` bits for each constraint. But this shouldn't cause any problems, as the 0 would have been caught by the prior check (`[x, 1, 0]`).
"""

# ╔═╡ 4a22d0ff-b089-4219-8117-598275e9a4cf
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

# ╔═╡ 17266853-7b9e-4836-8047-3ededf360bc5
#=╠═╡
function AbstractControllerAutomatonM(constraints::Vector{WeaklyHardConstraint})
	maxwindowsize = maximum(c -> c.windowsize, constraints)

	L = 2^maxwindowsize
	T = function (l, σ)
		l_new = (l << 1) & (L - 1) | σ
		l == 0 || !check_constraint(l_new, constraints) ? 0 : l_new
	end

	Automaton(L, 2, T, L-1, !iszero)
end
  ╠═╡ =#

# ╔═╡ b85e9591-8882-46ab-afef-8154a50e7281
let
	c1 = AbstractControllerAutomaton(WeaklyHardConstraint(1, 2))
	c2 = AbstractControllerAutomaton(WeaklyHardConstraint(1, 3))
end

# ╔═╡ 51b1941f-d4c7-44f0-9cc2-7cd69fa6a0f6
let
	c1 = AbstractControllerAutomaton(WeaklyHardConstraint(1, 2))
	[[l for l in 0:c1.L-1] [c1.T(l, σ) for l in 0:c1.L-1, σ in 0:c1.Σ-1]]
end

# ╔═╡ 97336b50-d02f-4150-8440-46b39caf80ce
let
	c2 = AbstractControllerAutomaton(WeaklyHardConstraint(1, 3))
	[[l for l in 0:c2.L-1] [c2.T(l, σ) for l in 0:c2.L-1, σ in 0:c2.Σ-1]]
end

# ╔═╡ a7eaa247-7557-40ac-b668-c53d3ab4dd73
md"""
Let's do some testing to make sure that the transition function of our `AbstractControllerAutomaton` is correct.
"""

# ╔═╡ f14406b4-cee3-48eb-b76e-27805b5073b3
# ╠═╡ disabled = true
#=╠═╡
let
	c3 = AbstractControllerAutomaton([
		WeaklyHardConstraint(1, 3)
		WeaklyHardConstraint(2, 2)
	])
	[[l for l in 0:c3.L-1] [c3.T(l, σ) for l in 0:c3.L-1, σ in 0:c3.Σ-1]]
end
  ╠═╡ =#

# ╔═╡ ca837ee8-8cdd-4d15-acd4-41d484c99d79
md"""
## Generating scheduler automaton

Before defining the generator function for the scheduler automaton we need to define a few functions.
"""

# ╔═╡ 9eb61af7-c961-4a7f-a6cc-ac70d99f19d1
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

# ╔═╡ 7fd78424-2c3f-49e1-a4ba-097762e4a48c
#=╠═╡
check_constraint(undigit([1, 0, 1]), [
	WeaklyHardConstraint(2, 3),
	WeaklyHardConstraint(1, 2),
	WeaklyHardConstraint(1, 1)
])
  ╠═╡ =#

# ╔═╡ 8a00a055-2f46-4788-bc0b-11c9b562f829
undigit([1, 0, 0])

# ╔═╡ 1fcf7a4b-a507-4bc1-909b-613849e336b2
"""
	digits_b2r(x[, pad])

A shortcut for `digits(x, base=2, pad=pad) |> reverse`
"""
function digits_b2r(x::Int, pad::Int=0)
	digits(x, base=2, pad=pad) |> reverse
end

# ╔═╡ c152607a-9f1d-446e-8d3d-a9e410255413
digits_b2r(5, 4)

# ╔═╡ 7d432d8e-7382-436b-a2f7-d98a637374c8
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

# ╔═╡ 7e719add-a193-4829-9c70-aceaabd0326a
state_separation(6, [1, 2])

# ╔═╡ 13863ec3-1805-406a-ac43-1ae3e43bfaed
md"""
And finally, the generator function.
"""

# ╔═╡ e00f079c-4be0-4876-bd59-4e760a19137b
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

# ╔═╡ 065bfbf2-743f-4e67-81a4-78a7ea2879de
SchedulerAutomaton(controllers::Vector{Automaton}; cores=1) =
	SchedulerAutomaton(controllers..., cores=cores)

# ╔═╡ 9bc2d9da-3c78-4d6e-8e7e-1f50e8824cc7
let
	controllers = AbstractControllerAutomaton.([
		WeaklyHardConstraint(1,7),
		WeaklyHardConstraint(3,7),
		WeaklyHardConstraint(2,6),
		WeaklyHardConstraint(1,4),
		WeaklyHardConstraint(2,5),
		WeaklyHardConstraint(4,5),
		WeaklyHardConstraint(3,6)
	])
	@time SchedulerAutomaton(controllers, cores=2)
end

# ╔═╡ 41f3a51b-4abe-4c3b-980c-80b46e265967
md"""
We test our implementation here (*to show the full result, run this notebook file in terminal and then inspect test1-3*)
"""

# ╔═╡ a7448e50-675c-4c6e-ad0c-5c1bf178e831
function scheduler_test(constraints::WeaklyHardConstraint...; cores=1)
	controllers = AbstractControllerAutomaton.(constraints)
	sa = SchedulerAutomaton(controllers..., cores=cores)
	headers = [missing, digits_b2r.(sa.Σ, sa.N)...]'
	current_states = [l for l in 0:sa.L-1]
	new_states = [sa.T(l, σ) for l in 0:sa.L-1, σ in sa.Σ]
	vcat(headers, map(l -> state_separation(l, sa.B), [current_states new_states]))
end

# ╔═╡ ea617897-7785-4b61-88d9-5b418de34332
test1 = let
	c1 = WeaklyHardConstraint(1, 2)
	c2 = WeaklyHardConstraint(2, 3)
	scheduler_test(c1, c2, cores=1)
end

# ╔═╡ dc40c20b-963b-4781-8dc5-192a893e532f
test2 = let
	c1 = WeaklyHardConstraint(1, 2)
	c2 = WeaklyHardConstraint(1, 3)
	scheduler_test(c1, c2, cores=2)
end;

# ╔═╡ 269d9339-8cbd-4d94-8d5c-c78c68c070d2
test3 = let
	c1 = WeaklyHardConstraint(1, 2)
	c2 = WeaklyHardConstraint(1, 2)
	c3 = WeaklyHardConstraint(2, 2)
	scheduler_test(c1, c2, c3, cores=1)
end;

# ╔═╡ b13fa861-9bce-4423-81e1-ceca425e1f2f
md"""
## Emptiness checking for scheduler automaton

Since the scheduler automaton is constructed so that a location is final if and only if all consisting locations for individual automata are final, we can check if a safe schedule exists by check for emptiness over a given number of steps, e.g., 100 steps. 

Furthermore, if there is a cycle detected within the automaton, then it is automatically considered success since the cycle can be followed in indefinitely many steps.
"""

# ╔═╡ be27305f-2dd2-4f0e-bf7a-a3d7589e3543
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

# ╔═╡ 2d1a79f8-b18e-41e5-9cb8-2a7327ff398b
let
	c1 = WeaklyHardConstraint(1, 2)
	c2 = WeaklyHardConstraint(1, 3)
	c3 = WeaklyHardConstraint(2, 3)
	path = isschedulable(c1, c2, c3, cores=2, steps=10)
	steps = map(l -> state_separation(l, [2, 3, 3], indigits=true), path) |> collect
end

# ╔═╡ 16919960-6191-47bc-a675-f1b0e78ac077
md"""
## Examples from real controllers

In this section we import the result from running the bounded tree deviation estimation algorithm on a few real world controllers and investigate the schedulability of them.
"""

# ╔═╡ 3e264d9b-606e-48be-9989-0d8659bc8e05
"""
	load_result(path, threashold)

Load the result (in `.csv` format) from the file system, and present the result according to the threashold value.
"""
function load_result(path, threshold)
	d = readdlm(path, ',', Float64)
	d = round.(d; digits=2)
	w = size(d, 2)
	deviations = d[1:w,:]
	indices    = d[w+1:2w,:]

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
	constraints
end

# ╔═╡ 60c97203-0b91-4677-9d98-8f5f9a861aaf
@bind threshold_1_2 Slider(1:0.01:2, default=1.8, show_value=true)

# ╔═╡ f077e353-2271-4fea-8b78-b4db76854b74
rc = load_result("../experiments/data/closeperiod/0.02s_0.02s/RC_HoldAndKill_n5_t100.csv", threshold_1_2)

# ╔═╡ c773d228-5bd4-4fc9-903f-2df96f0827b8
@bind threshold Slider(0:0.1:10, default=10, show_value=true)

# ╔═╡ 11963fb9-998b-4cc9-9f1d-c30f273ba906
f1 = load_result("../experiments/data/closeperiod/0.02s_0.02s/F1_HoldAndKill_n5_t100.csv", 5)

# ╔═╡ 2e2a92c3-442c-4061-a502-79357805139d
@bind threshold_500 Slider(0:100, default=20, show_value=true)

# ╔═╡ ee91e037-e94b-48bf-aef3-7dfe02f2bf3e
dcm = load_result("../experiments/data/closeperiod/0.02s_0.02s/DCM_HoldAndKill_n5_t100.csv", 0.24)

# ╔═╡ 5573a92e-224c-4850-9dea-2e122195e4de
css = load_result("../experiments/data/closeperiod/0.02s_0.02s/CSS_HoldAndKill_n5_t100.csv", 2.37)

# ╔═╡ 4ec2ad2d-4a3b-472f-93f6-26a913453e59
cc2 = load_result("../experiments/data/closeperiod/0.02s_0.02s/CC2_HoldAndKill_n5_t100.csv", 5)

# ╔═╡ 9b74b2be-72d1-4ed4-ba75-6731361f8fe4
# urgency+largest window counter example = [(1, 4), (3, 6), (1, 5)]
let 
	c3 = WeaklyHardConstraint(1, 4)
	c4 = WeaklyHardConstraint(3, 6)
	c5 = WeaklyHardConstraint(1, 5)

	path = isschedulable(c3, c4, c5, cores=1)
	if length(path) > 0
		@info "Found results" c3 c4 c5 map(l -> state_separation(l, [c.windowsize for c in [c3, c4, c5]], indigits=true), path) |> collect
	end
end

# ╔═╡ 701e5540-120d-4f8c-a693-bb1ac56f40ed
# for c1 in ap, c2 in es, c3 in f1
# 	path = isschedulable(c1, c2, c3, cores=2)
# 	if length(path) > 0
# 		@info "Found results" c1 c2 c3 map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3]], indigits=true), reverse(path)) |> collect
# 	end
# end
for c1 in rc, c2 in f1, c3 in dcm, c4 in css, c5 in cc2
	path = isschedulable(c1, c2, c3, c4, c5, cores=2)
	if length(path) > 0
		@info "Found results" c1 c2 c3 c4 c5 map(l -> state_separation(l, [c.windowsize for c in [c1, c2, c3, c4, c5]], indigits=true), path) |> collect
		break
	end
end

# in figure show deviation graph

# ╔═╡ 84add455-d94f-49af-8477-bd60965e52d1
# ╠═╡ disabled = true
#=╠═╡
isschedulable(
	WeaklyHardConstraint(1, 7),
	WeaklyHardConstraint(2, 7),
	WeaklyHardConstraint(3, 7),
	WeaklyHardConstraint(4, 7),
	cores=2
)
  ╠═╡ =#

# ╔═╡ 014a97e6-4f86-4c58-9a55-84734f01a9ab
# ╠═╡ disabled = true
#=╠═╡
let
	allconstraints = [WeaklyHardConstraint(h, w) for w=2:5 for h=1:w-1]
	for numcontrollers=2:4
		combo = Iterators.product(repeat([allconstraints], numcontrollers)...)
		for cores=1:numcontrollers-1
			for constraints in combo
				isschedulable(constraints..., cores=cores)
			end
		end
		@info "Finished numcontrollers=$(numcontrollers)"
	end
end
  ╠═╡ =#

# ╔═╡ f830a585-647c-4adf-910e-18f870ff73cf
# ╠═╡ disabled = true
#=╠═╡
let
	allconstraints = [WeaklyHardConstraint(h, w) for w=6 for h=1:w-1]
	for numcontrollers=2:4
		combo = Iterators.product(repeat([allconstraints], numcontrollers)...)
		for cores=1:numcontrollers-1
			for constraints in combo
				isschedulable(constraints..., cores=cores)
			end
		end
		@info "Finished numcontrollers=$(numcontrollers)"
	end
end
  ╠═╡ =#

# ╔═╡ 75c5f9ce-927e-4944-bcdc-38e72abbd982
Threads.nthreads()

# ╔═╡ ee68db72-8387-4dfa-bf1f-630a3b44170d
# ╠═╡ disabled = true
#=╠═╡
let
	allconstraints = [WeaklyHardConstraint(h, w) for w=2:5 for h=1:w-1]
	for numcontrollers=5
		combo = Iterators.product(repeat([allconstraints], numcontrollers)...)
		for cores=1:numcontrollers-1
			for constraints in combo
				isschedulable(constraints..., cores=cores)
			end
		end
		@info "Finished numcontrollers=$(numcontrollers)"
	end
end
  ╠═╡ =#

# ╔═╡ 128a2e0d-0c90-4eb6-ad2b-495717f7d85d
# ╠═╡ disabled = true
#=╠═╡
let
	allconstraints = [WeaklyHardConstraint(h, w) for w=6 for h=1:w-1]
	for numcontrollers=5
		combo = Iterators.product(repeat([allconstraints], numcontrollers)...)
		for cores=1:numcontrollers-1
			for constraints in combo
				isschedulable(constraints..., cores=cores)
			end
		end
		@info "Finished numcontrollers=$(numcontrollers)"
	end
end
  ╠═╡ =#

# ╔═╡ a46320cb-7297-443d-b16c-26f6dc2efb43
let
	allconstraints = [WeaklyHardConstraint(h, w) for w=2:6 for h=1:w-1]
	numcontrollers = 3
	size(
		Iterators.product(
			repeat([allconstraints], numcontrollers)...
		) |> collect
	)
end

# ╔═╡ Cell order:
# ╟─191d8766-d7d4-11ec-2c79-c5eee06d2d64
# ╠═3fbc35e1-e930-4a60-9eb8-283d733901d3
# ╠═b6f14d07-77aa-4948-af8e-7d6062d6c38f
# ╟─df7c19a5-69a4-4bd6-a157-572829150052
# ╠═2c3b7ca2-cda4-4274-a14c-8a57c16d15cd
# ╠═3af3641c-687b-42bb-bfea-0f73425cb827
# ╟─f1de5c80-c5bf-4d8d-bbb1-78ce7d675ff5
# ╠═e7bf324a-c672-4161-84c9-4993efd8a0af
# ╟─b39b4969-03df-41d9-99b6-fb95d98a08ba
# ╠═ce56b597-02c5-4d3d-9712-762b493ae7a4
# ╠═d8f5a9cf-d535-40c4-b193-54d1b8a00e50
# ╟─5da180a4-7a0e-463d-bbb3-74e83a23eb47
# ╠═7fd78424-2c3f-49e1-a4ba-097762e4a48c
# ╠═4a22d0ff-b089-4219-8117-598275e9a4cf
# ╠═17266853-7b9e-4836-8047-3ededf360bc5
# ╠═b85e9591-8882-46ab-afef-8154a50e7281
# ╠═51b1941f-d4c7-44f0-9cc2-7cd69fa6a0f6
# ╠═97336b50-d02f-4150-8440-46b39caf80ce
# ╟─a7eaa247-7557-40ac-b668-c53d3ab4dd73
# ╠═f14406b4-cee3-48eb-b76e-27805b5073b3
# ╟─ca837ee8-8cdd-4d15-acd4-41d484c99d79
# ╠═9eb61af7-c961-4a7f-a6cc-ac70d99f19d1
# ╠═8a00a055-2f46-4788-bc0b-11c9b562f829
# ╠═1fcf7a4b-a507-4bc1-909b-613849e336b2
# ╠═c152607a-9f1d-446e-8d3d-a9e410255413
# ╠═7d432d8e-7382-436b-a2f7-d98a637374c8
# ╠═7e719add-a193-4829-9c70-aceaabd0326a
# ╟─13863ec3-1805-406a-ac43-1ae3e43bfaed
# ╠═e00f079c-4be0-4876-bd59-4e760a19137b
# ╠═065bfbf2-743f-4e67-81a4-78a7ea2879de
# ╠═9bc2d9da-3c78-4d6e-8e7e-1f50e8824cc7
# ╟─41f3a51b-4abe-4c3b-980c-80b46e265967
# ╠═a7448e50-675c-4c6e-ad0c-5c1bf178e831
# ╠═ea617897-7785-4b61-88d9-5b418de34332
# ╠═dc40c20b-963b-4781-8dc5-192a893e532f
# ╠═269d9339-8cbd-4d94-8d5c-c78c68c070d2
# ╟─b13fa861-9bce-4423-81e1-ceca425e1f2f
# ╠═be27305f-2dd2-4f0e-bf7a-a3d7589e3543
# ╠═2d1a79f8-b18e-41e5-9cb8-2a7327ff398b
# ╟─16919960-6191-47bc-a675-f1b0e78ac077
# ╠═b1c60c0a-e804-461b-a026-f173d963599c
# ╠═3e264d9b-606e-48be-9989-0d8659bc8e05
# ╠═60c97203-0b91-4677-9d98-8f5f9a861aaf
# ╠═f077e353-2271-4fea-8b78-b4db76854b74
# ╠═c773d228-5bd4-4fc9-903f-2df96f0827b8
# ╠═11963fb9-998b-4cc9-9f1d-c30f273ba906
# ╠═2e2a92c3-442c-4061-a502-79357805139d
# ╠═ee91e037-e94b-48bf-aef3-7dfe02f2bf3e
# ╠═5573a92e-224c-4850-9dea-2e122195e4de
# ╠═4ec2ad2d-4a3b-472f-93f6-26a913453e59
# ╠═9b74b2be-72d1-4ed4-ba75-6731361f8fe4
# ╠═701e5540-120d-4f8c-a693-bb1ac56f40ed
# ╠═84add455-d94f-49af-8477-bd60965e52d1
# ╠═014a97e6-4f86-4c58-9a55-84734f01a9ab
# ╠═f830a585-647c-4adf-910e-18f870ff73cf
# ╠═75c5f9ce-927e-4944-bcdc-38e72abbd982
# ╠═ee68db72-8387-4dfa-bf1f-630a3b44170d
# ╠═128a2e0d-0c90-4eb6-ad2b-495717f7d85d
# ╠═a46320cb-7297-443d-b16c-26f6dc2efb43
