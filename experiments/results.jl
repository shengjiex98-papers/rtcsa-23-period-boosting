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

# ╔═╡ 41a85aac-d357-11ec-3d01-71de103a66ec
begin
	using JLD
	using DelimitedFiles
	using PlutoUI
	using Logging
end

# ╔═╡ d6708b2d-0ccd-472b-b994-a175b5960d29
function getdev(path; sigdigits=2)
	d = readdlm(path, ',', Float64)
	d = round.(d; sigdigits=sigdigits)
	w = size(d, 2)
	dev = d[1:w,:]
	ind = d[w+1:2w,:]

	dev, ind, w
end

# ╔═╡ e2a0c124-bf14-4fbc-bddb-959d1220581b
# Logging.loglevel(Logging.Warn)
function summer(path, threshold_percent; sigdigits=2)
	deviations, indices, w = getdev(path, sigdigits=sigdigits)
	
	min_v = minimum(x -> isnan(x) ? Inf  : x, deviations)
	max_v = maximum(x -> isnan(x) ? -Inf : x, deviations)
	threshold = min_v + (max_v - min_v) * threshold_percent

	to_show = fill(NaN, size(deviations))

	min_hits = 1
	for window_size in 2:w+1
		while min_hits < window_size
			if deviations[window_size-1, min_hits] <= threshold
				to_show[window_size-1, min_hits] = deviations[window_size-1, min_hits]
				break
			end
			min_hits += 1
		end
	end
	
	# @info deviations
	# @info deviations indices
	display("text/plain", [deviations; indices])
	to_show
end

# ╔═╡ acb69bf3-7f44-4e49-9e04-619e455db184
function compare_gain(path1, path2, threshold_percent)
	dev1, ind1, w = getdev(path1, sigdigits=3)
	dev2, ind2, w = getdev(path2, sigdigits=3)

	min_v = minimum(x -> isnan(x) ? Inf  : x, vcat(dev1, dev2))
	max_v = maximum(x -> isnan(x) ? -Inf : x, vcat(dev1, dev2))
	threshold = min_v + (max_v - min_v) * threshold_percent
	
	to_show = fill(' ', (w, w))

	min_hits = 1
	for window_size in 2:w+1
		for min_hits in 1:window_size-1
			d1 = dev1[window_size-1, min_hits]
			d2 = dev2[window_size-1, min_hits]
			if d1 <= threshold && d2 <= threshold
				sy = 's'
			elseif d1 > threshold && d2 > threshold
				sy = 'x'
			elseif d1 <= threshold && d2 > threshold
				sy = '1'
			elseif d1 > threshold && d2 <= threshold
				sy = '2'
			else
				sy = '?'
			end
			to_show[window_size-1, min_hits] = sy
			# min_hits += 1
		end
	end
	
	# @info deviations
	# @info deviations indices
	# display("text/plain", [dev1; ind1])
	@info threshold
	to_show
end

# ╔═╡ 5f80af0d-d220-4f22-8507-c8dae4335433
@bind RC Slider(0:0.01:1, show_value=true)

# ╔═╡ 6bb79ca7-2dc1-4ce8-aceb-621356a96bca
let
	path1 = "./data/finalize/0.015s_0.015s/DCM_100_n10_t100.csv"
	path2 = "./data/finalize/0.015s_0.023s/DCM_100_n10_t100.csv"
	compare_gain(path1, path2, RC)
end

# ╔═╡ 66731c94-e83e-414c-80eb-b37d8ef69724
md"""
## 20 ms period w/ 20 ms discretization
"""

# ╔═╡ 981a2d06-32a8-4fe8-8b10-c9525f72ddf3
@bind perc Slider(0:0.01:1, show_value=true)

# ╔═╡ 1406c32e-ce99-4de6-aa29-1d8ee2cf7744
# ╠═╡ disabled = true
#=╠═╡
let
	summer("data/closeperiod/0.018s_0.02s/RC_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	summer("data/closeperiod/0.018s_0.02s/DCM_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	summer("data/closeperiod/0.018s_0.02s/CSS_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	# summer("data/closeperiod/0.018s_0.02s/CC1_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	summer("data/closeperiod/0.018s_0.02s/CC2_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	md"""
	## 18 ms period w/ 20 ms discretization
	"""
end
  ╠═╡ =#

# ╔═╡ 9b273a22-552e-4b35-974e-7675aebf5fa4
let
	summer("data/closeperiod/0.015s_0.02s/RC_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	summer("data/closeperiod/0.015s_0.02s/DCM_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	summer("data/closeperiod/0.015s_0.02s/CSS_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	# summer("data/closeperiod/0.015s_0.02s/CC1_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	summer("data/closeperiod/0.015s_0.02s/CC2_HoldAndKill_n5_t100.csv", perc, sigdigits=3)
	md"""
	## 15 ms period w/ 20 ms discretization
	"""
end

# ╔═╡ b73c20d7-6271-4533-b986-54289e9b357b
summer("data/common_period/0.02s_0.02s/RC_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ c4ab5843-067c-44a6-aec9-7f0290f4853a
summer("data/common_period/0.02s_0.02s/DCM_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 485404b8-ad5c-4e52-9806-dbc2b6623cef
summer("data/common_period/0.02s_0.02s/CSS_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 992a17e8-4cc5-4cdb-8ac2-8489d16ad205
summer("data/common_period/0.02s_0.02s/CC1_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 3c96a791-8722-4a3f-b071-6231cfb95aa9
summer("data/default/0.02s_0.02s/CC2_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 0a47ccae-a40d-411b-bb26-3953581997af
md"""
## 23 ms period w/ 20 ms discretization
"""

# ╔═╡ ca8f8a4c-ad2b-40a2-8e39-481080d07e44
summer("data/closeperiod/0.023s_0.02s/DCM_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 7b0c1274-4208-4c1f-b12d-b5647b50f98f
summer("data/closeperiod/0.023s_0.02s/CSS_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 6124477b-971e-4bdc-9210-1af57b485877
summer("data/closeperiod/0.023s_0.02s/CC1_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ a52aefd5-87d7-490e-9445-3bb230e5d684
summer("data/closeperiod/0.023s_0.02s/CC2_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ d4904ddf-4d5e-426a-80c1-71152ebf4f3b
md"""
## 25 ms period w/ 20 ms discretization
"""

# ╔═╡ 96d1e07a-d0e8-4909-a56f-2aaf273e8e7c
summer("data/closeperiod/0.025s_0.02s/DCM_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 5e4d593c-1dcf-4f11-8762-6726d35144c4
summer("data/closeperiod/0.025s_0.02s/CSS_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 87634f75-0eb8-4fd3-b69e-e9d94de4d80f
summer("data/closeperiod/0.025s_0.02s/CC1_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 8470c142-1539-4a2d-bf26-8698adc6e647
summer("data/closeperiod/0.025s_0.02s/CC2_HoldAndKill_n5_t100.csv", perc, sigdigits=3)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
JLD = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
Logging = "56ddb016-857b-54e1-b83d-db4d58db5568"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
JLD = "~0.13.1"
PlutoUI = "~0.7.38"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Blosc]]
deps = ["Blosc_jll"]
git-tree-sha1 = "310b77648d38c223d947ff3f50f511d08690b8d5"
uuid = "a74b3585-a348-5f62-a45c-50e91977d574"
version = "0.7.3"

[[Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "91d6baa911283650df649d0aea7c28639273ae7b"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.1+0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "94f5101b96d2d968ace56f7f2db19d0a5f592e28"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.15.0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[H5Zblosc]]
deps = ["Blosc", "HDF5"]
git-tree-sha1 = "26b22c9039b039e29ec4f4f989946de722e87ab9"
uuid = "c8ec2601-a99c-407f-b158-e79c03c2f5f7"
version = "0.1.0"

[[HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires"]
git-tree-sha1 = "9ffc57b9bb643bf3fce34f3daf9ff506ed2d8b7a"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.10"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "c003b31e2e818bc512b0ff99d7dce03b0c1359f5"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JLD]]
deps = ["Compat", "FileIO", "H5Zblosc", "HDF5", "Printf"]
git-tree-sha1 = "cd46c18390e9bbc37a2098dfb355ec5f18931900"
uuid = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
version = "0.13.2"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═41a85aac-d357-11ec-3d01-71de103a66ec
# ╠═d6708b2d-0ccd-472b-b994-a175b5960d29
# ╠═e2a0c124-bf14-4fbc-bddb-959d1220581b
# ╠═acb69bf3-7f44-4e49-9e04-619e455db184
# ╠═5f80af0d-d220-4f22-8507-c8dae4335433
# ╠═6bb79ca7-2dc1-4ce8-aceb-621356a96bca
# ╟─1406c32e-ce99-4de6-aa29-1d8ee2cf7744
# ╠═9b273a22-552e-4b35-974e-7675aebf5fa4
# ╠═66731c94-e83e-414c-80eb-b37d8ef69724
# ╠═981a2d06-32a8-4fe8-8b10-c9525f72ddf3
# ╠═b73c20d7-6271-4533-b986-54289e9b357b
# ╠═c4ab5843-067c-44a6-aec9-7f0290f4853a
# ╠═485404b8-ad5c-4e52-9806-dbc2b6623cef
# ╠═992a17e8-4cc5-4cdb-8ac2-8489d16ad205
# ╠═3c96a791-8722-4a3f-b071-6231cfb95aa9
# ╟─0a47ccae-a40d-411b-bb26-3953581997af
# ╠═ca8f8a4c-ad2b-40a2-8e39-481080d07e44
# ╠═7b0c1274-4208-4c1f-b12d-b5647b50f98f
# ╠═6124477b-971e-4bdc-9210-1af57b485877
# ╠═a52aefd5-87d7-490e-9445-3bb230e5d684
# ╟─d4904ddf-4d5e-426a-80c1-71152ebf4f3b
# ╠═96d1e07a-d0e8-4909-a56f-2aaf273e8e7c
# ╠═5e4d593c-1dcf-4f11-8762-6726d35144c4
# ╠═87634f75-0eb8-4fd3-b69e-e9d94de4d80f
# ╠═8470c142-1539-4a2d-bf26-8698adc6e647
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
