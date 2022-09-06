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
	using PlutoUI
	using DelimitedFiles
	using DataStructures
end

# ╔═╡ d2fe79ad-3b5f-461c-be6f-f3e7aaafe8fa
TableOfContents()

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
	function load_result(path, threshold_ratio)
		d = readdlm(path, ',', Float64)
		d = round.(d; sigdigits=3)
		w = size(d, 2)
		deviations = d[1:w,:]
		indices    = d[w+1:2w,:]
		# @info deviations

		min_v = minimum(x -> isnan(x) ? Inf  : x, deviations)
		max_v = maximum(x -> isnan(x) ? -Inf : x, deviations)
		threshold = round(min_v + (max_v - min_v) * threshold_ratio, sigdigits=3)
	
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
| RC     | 10ms | 23ms   | 0.0159        |
| F1     | 13ms | 20ms   | 1520.0        |
| DCM    | 12ms | 23ms   | 0.0336        |
| CSS    | 10ms | 27ms   | 57.5          |
| CC     | 15ms | 28ms   | 1.05          |
Utilization $$U \approx 2.51 > 1$$
"""

# ╔═╡ 470a7307-82a4-4a30-a308-6d8ade16a2f4
prefix = "../experiments/data/common_period/"

# ╔═╡ 57c2fd32-48b2-4a53-9bb3-96f28ef6c70b
@bind threshold_rc Slider(0:0.01:1, default=0.5, show_value=true)

# ╔═╡ f17e5f30-8dbb-4079-a22b-d79b05aedbed
rc = load_result(prefix * "0.028s_0.023s/RC_HoldAndKill_n5_t100.csv", threshold_rc)

# ╔═╡ 421842ca-442c-499f-8ff5-1ac46d53e67e
@bind threshold_f1 Slider(0:0.01:1, default=0.5, show_value=true)

# ╔═╡ e3f33864-27e5-4bff-be7f-942d3a8cdaf8
f1 = load_result("../experiments/data/common_period/0.028s_0.02s/F1_HoldAndKill_n5_t100.csv", threshold_f1)

# ╔═╡ f37293d5-3ba0-4abf-898a-8d500d11302e
@bind threshold_dcm Slider(0:0.01:1, default=0.5, show_value=true)

# ╔═╡ 4d862328-ca2f-4656-a296-e08237fd3150
dcm = load_result(prefix * "0.028s_0.023s/DCM_HoldAndKill_n5_t100.csv", threshold_dcm)

# ╔═╡ b3dd90f5-c87d-4d67-8094-3eaf836033f0
@bind threshold_css Slider(0:0.01:1, default=0.5, show_value=true)

# ╔═╡ 076dc76e-b93c-4642-a189-5bf233ba4eb0
css = load_result(prefix * "0.028s_0.027s/CC2_HoldAndKill_n5_t100.csv", threshold_css)

# ╔═╡ adb5ed1f-d1ab-4eda-a362-5fd26e98b6c7
@bind threshold_cc2 Slider(0:0.01:1, default=0.5, show_value=true)

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DataStructures = "~0.18.13"
PlutoUI = "~0.7.39"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0"
manifest_format = "2.0"
project_hash = "2ff104bdee91d2b2e23bdc6364e673ea35697286"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "5856d3031cdb1f3b2b6340dfdc66b6d9a149a374"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.2.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

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
# ╟─421842ca-442c-499f-8ff5-1ac46d53e67e
# ╠═e3f33864-27e5-4bff-be7f-942d3a8cdaf8
# ╟─f37293d5-3ba0-4abf-898a-8d500d11302e
# ╠═4d862328-ca2f-4656-a296-e08237fd3150
# ╟─b3dd90f5-c87d-4d67-8094-3eaf836033f0
# ╠═076dc76e-b93c-4642-a189-5bf233ba4eb0
# ╠═adb5ed1f-d1ab-4eda-a362-5fd26e98b6c7
# ╠═71875c91-2167-4fa6-a26d-91f779269796
# ╠═29e7a702-5ef4-4f1f-bfcb-5a0b69621749
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
