### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ ba91b96e-1c2a-11ee-35f8-2183fe15f018
begin
	import Pkg
	Pkg.activate("..")

	using RealTimeScheduling
	using ControlTimingSafety
	using ControlSafetyBench
	
	using Plots, LaTeXStrings, PlutoUI
	using DelimitedFiles

	push!(LOAD_PATH, "../lib")
	# include("../lib/Experiments.jl")
	using Experiments
	
	TableOfContents()
end

# ╔═╡ 6e8ceb64-c5a7-4dee-8078-ce5aed37a169
md"""
# Schedule Synthesis Analysis
"""

# ╔═╡ 5470e671-b29a-44e1-9bdd-e2cb20937c52
md"""
## Reading Data from Files
"""

# ╔═╡ 4eb354c4-459c-402b-9bcd-2c4bd00518bd
systems = [:RC, :F1, :DC, :CS, :CC]

# ╔═╡ 97f1fd06-cebf-427d-a107-75edf27cc544
ps = [15, 28, 40]

# ╔═╡ 02ab2b94-f667-483c-9d20-4e12f836436d
rs = ["", "-RE"]

# ╔═╡ c0fb0533-9826-4d78-8ae2-4ef811fa3630
data = let
	dir = "../data/"
	rdev(sys, p, r) = readdlm("$dir$sys-$p$r.tsv", '\t', Float64, '\n')
	list = map(Iterators.product(systems, ps, rs)) do (sys, p, r)
		Symbol(sys, p, r) => rdev(sys, p, r)
	end
	(; list...)
end

# ╔═╡ f25c2946-3a77-4c28-b720-afbe11047721
md"""
## Finding safe constraints
"""

# ╔═╡ e96a3c24-e916-4b60-b380-59ea113f9810
saf = Dict()

# ╔═╡ f4183950-4f18-420b-b592-7669e3343194
@bind tRC Slider(1:0.01:10, default=5.8, show_value=true)

# ╔═╡ b65da002-1de4-485e-93a7-06722af981c3
@bind tF1 Slider(0.1:0.001:1, default=0.37, show_value=true)

# ╔═╡ fce45e39-6f86-4ef5-a0d2-15b9c2bdd434
@bind tDC Slider(0.001:0.001:0.5, default=0.18, show_value=true)

# ╔═╡ 197931fb-d60a-4434-8d77-3f079793ff88
@bind tCS Slider(0.1:0.005:5, default=2.1, show_value=true)

# ╔═╡ 37087901-56b9-4e2b-beef-3ebf5d0c85a8
@bind tCC Slider(0.05:0.005:1, default=0.48, show_value=true)

# ╔═╡ e32ad5df-0692-4870-8ff8-e3a4c106060d
function find_safe!(saf, data, sys::Symbol, t::Real)
	maxwindow = size(data[Symbol(sys, 15)], 1)
	res = fill("", maxwindow, maxwindow)
	for p in ps, r in rs
		id = Symbol(sys, p, r)
		saf[id] = Vector{MeetAny{Int64}}()
		for window = 1:maxwindow, meet = 1:window
			if data[id][window, meet] < t 
				res[window, meet] *= "$(p÷10)"
				push!(saf[id], MeetAny(meet, window))
			else
				res[window, meet] *= " "
			end
		end
	end
	@info res
	saf
end

# ╔═╡ c3d7f90c-6962-4eca-bbc0-20d0af4125aa
r1 = find_safe!(saf, data, :RC, tRC);

# ╔═╡ 045260ed-6a24-4359-ae7c-ff12dcc8a20c
r2 = find_safe!(saf, data, :F1, tF1);

# ╔═╡ e2493314-841f-42d1-9a3d-7e334650e9ec
r3 = find_safe!(saf, data, :DC, tDC);

# ╔═╡ 7f5b7f07-e5fe-4425-a59d-3b37bcc53778
r4 = find_safe!(saf, data, :CS, tCS);

# ╔═╡ a402a64a-2a42-42a2-9d49-b4a225119865
r5 = find_safe!(saf, data, :CC, tCC);

# ╔═╡ 0154f2df-6e1a-4d66-b3a9-08b96d682608
md"""
Now, use the safe constraints to attempt schedule synthesis
"""

# ╔═╡ 28d67e5c-378a-4863-8cab-731cbea17f63
s15 = let
	r1, r2, r3, r4, r5
	consn = map(systems) do sys
		saf[Symbol(sys, 15)]
	end
	consr = map(systems) do sys
		saf[Symbol(sys, 15, "-RE")]
	end
	schedule_xghtc(consn, 100, slotsize=1),
	schedule_xghtc(consr, 100, slotsize=1)
end

# ╔═╡ 9c76bed6-0e5a-4219-8881-b22d418b05a2
s28 = let
	r1, r2, r3, r4, r5
	consn = map(systems) do sys
		saf[Symbol(sys, 28)]
	end
	consr = map(systems) do sys
		saf[Symbol(sys, 28, "-RE")]
	end
	schedule_xghtc(consn, 100, slotsize=2, fullpath=false),
	schedule_xghtc(consr, 100, slotsize=2, fullpath=false)
end

# ╔═╡ 37da5df9-103b-413d-9de0-9693f1eb7263
s40 = let
	r1, r2, r3, r4, r5
	consn = map(systems) do sys
		saf[Symbol(sys, 40)]
	end
	consr = map(systems) do sys
		saf[Symbol(sys, 40, "-RE")]
	end
	schedule_xghtc(consn, 100, slotsize=3, fullpath=true),
	schedule_xghtc(consr, 100, slotsize=3, fullpath=true)
end

# ╔═╡ 516e761e-d193-45a5-9931-850d99c6f5fa
md"""
## Deviation Visualization
"""

# ╔═╡ 4858b2e4-c2dd-4ab1-bfbb-862d16802bf5
md"""
### Camera Settings

|   |   |
|--:|:--|
| Azimuth 	| $(@bind az Slider(-180:1:180, default=35, show_value=true)) |
| Elevation | $(@bind el Slider(-90:1:90,   default=25, show_value=true)) |
"""

# ╔═╡ Cell order:
# ╟─6e8ceb64-c5a7-4dee-8078-ce5aed37a169
# ╠═ba91b96e-1c2a-11ee-35f8-2183fe15f018
# ╟─5470e671-b29a-44e1-9bdd-e2cb20937c52
# ╠═4eb354c4-459c-402b-9bcd-2c4bd00518bd
# ╠═97f1fd06-cebf-427d-a107-75edf27cc544
# ╠═02ab2b94-f667-483c-9d20-4e12f836436d
# ╠═c0fb0533-9826-4d78-8ae2-4ef811fa3630
# ╟─f25c2946-3a77-4c28-b720-afbe11047721
# ╠═e96a3c24-e916-4b60-b380-59ea113f9810
# ╠═f4183950-4f18-420b-b592-7669e3343194
# ╠═c3d7f90c-6962-4eca-bbc0-20d0af4125aa
# ╠═b65da002-1de4-485e-93a7-06722af981c3
# ╠═045260ed-6a24-4359-ae7c-ff12dcc8a20c
# ╠═fce45e39-6f86-4ef5-a0d2-15b9c2bdd434
# ╠═e2493314-841f-42d1-9a3d-7e334650e9ec
# ╠═197931fb-d60a-4434-8d77-3f079793ff88
# ╠═7f5b7f07-e5fe-4425-a59d-3b37bcc53778
# ╠═37087901-56b9-4e2b-beef-3ebf5d0c85a8
# ╠═a402a64a-2a42-42a2-9d49-b4a225119865
# ╠═e32ad5df-0692-4870-8ff8-e3a4c106060d
# ╟─0154f2df-6e1a-4d66-b3a9-08b96d682608
# ╠═28d67e5c-378a-4863-8cab-731cbea17f63
# ╠═9c76bed6-0e5a-4219-8881-b22d418b05a2
# ╠═37da5df9-103b-413d-9de0-9693f1eb7263
# ╟─516e761e-d193-45a5-9931-850d99c6f5fa
# ╠═4858b2e4-c2dd-4ab1-bfbb-862d16802bf5
