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

# ╔═╡ 41a85aac-d357-11ec-3d01-71de103a66ec
begin
	using JLD
	using DelimitedFiles
	using PlutoUI
end

# ╔═╡ 981a2d06-32a8-4fe8-8b10-c9525f72ddf3
@bind threshold_1_2 Slider(1:0.01:2, default=1.8, show_value=true)

# ╔═╡ 4e207296-a683-4f9c-863b-3205fdda2fbd
@bind threshold_0_01 Slider(0:0.001:0.1, default=0.05, show_value=true)

# ╔═╡ 9a5471c8-1d00-48ff-a74f-67448835a50c
@bind threshold Slider(0:0.1:15, default=10, show_value=true)

# ╔═╡ b9e019af-d155-4f1c-9de3-712f1fd02011
@bind threshold_ES Slider(0:0.0001:0.3, default=0.05, show_value=true)

# ╔═╡ d5d11eb5-8f05-4190-aacb-dca923fa894c
@bind threshold_F1 Slider(0:0.01:10, default=1, show_value=true)

# ╔═╡ 203ebc89-e958-441f-8788-75d1bc945a3f
@bind threshold_500 Slider(0:500, default=20, show_value=true)

# ╔═╡ db475d07-ed0e-478c-8fb3-5975719266ab
@bind threshold_AP Slider(0:0.1:5, default=4, show_value=true)

# ╔═╡ a9578595-f100-44dc-89b0-a8926c979acb
function summer(path, threshold; sigdigits=2)
	d = readdlm(path, ',', Float64)
	d = round.(d; sigdigits=sigdigits)
	w = size(d, 2)
	deviations = d[1:w,:]
	indices    = d[w+1:2w,:]
	
	# min_v = minimum(x -> isnan(x) ? Inf  : x, deviations)
	# max_v = maximum(x -> isnan(x) ? -Inf : x, deviations)
	# slider = Slider(0:1000, default=100, show_value=true)
	# @bind threshold slider

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
	
	@info deviations indices
	to_show
end

# ╔═╡ c4ab5843-067c-44a6-aec9-7f0290f4853a
summer("data/parallel/RC_HoldAndKill.csv", threshold_1_2, sigdigits=3)

# ╔═╡ 9665c251-7374-436f-b1fd-c8b5cdef8270
summer("data/period_0.02/RC_HoldAndKill.csv", threshold_0_01, sigdigits=3)

# ╔═╡ 37332697-7b8a-4f6c-a56b-6c29c9ea273c
summer("data/short/ES_HoldAndKill.csv", threshold, sigdigits=2)

# ╔═╡ 3fbd754c-7ce1-4d2c-83c0-e54fc451ed0e
summer("data/parallel/F1_HoldAndKill.csv", threshold, sigdigits=2)

# ╔═╡ 4fd31c79-c8fe-4d3e-b5a8-e111ae748277
summer("data/period_0.02/ES_HoldAndKill.csv", threshold_ES, sigdigits=3)

# ╔═╡ 13ec75c8-874e-4f08-8b8a-416d9f46613c
summer("data/period_0.02/F1_HoldAndKill.csv", threshold_F1, sigdigits=3)

# ╔═╡ 067d97f9-2599-4b64-8791-2d5d9fc49a07
summer("data/short/AP_HoldAndKill.csv", threshold_500)

# ╔═╡ f4c0af9a-57bb-4cdc-9f8d-ba6113e59f72
summer("data/period_0.02/AP_HoldAndKill.csv", threshold_AP)

# ╔═╡ 09ef9182-5350-4347-82e8-fc4cf13c3976
@bind t_dcm Slider(0:0.01:0.5, show_value=true)

# ╔═╡ 596ae945-2d62-4b2a-8349-6c5f416d1dda
summer("data/default/DCM_HoldAndKill_n5_t100.csv", t_dcm)

# ╔═╡ f1b3a897-70d5-4c68-98f7-6727e72fa442
@bind t_css Slider(0:0.01:12, show_value=true)

# ╔═╡ d8029862-0259-4982-94e8-496b4239063e
summer("data/default/CSS_HoldAndKill_n5_t100.csv", t_css)

# ╔═╡ ceb12397-a119-485a-a7e3-fdabfa4ae2f3
@bind t_ewb Slider(0:0.01:12, show_value=true)

# ╔═╡ 379c5629-e8df-4188-9716-b0f0db996119
summer("data/default/EWB_HoldAndKill_n5_t100.csv", t_ewb)

# ╔═╡ 692cea24-bf1a-4275-b388-8f2fc7ea610e
@bind t_cc1 Slider(0:0.01:12, show_value=true)

# ╔═╡ aba4fefb-f8db-4d37-a287-aa0e6352b4ff
summer("data/default/CC1_HoldAndKill_n5_t100.csv", t_cc1)

# ╔═╡ edfe4eb1-a007-47d3-a195-704479587365
@bind t_cc2 Slider(0:0.01:3, show_value=true)

# ╔═╡ 93d654ad-357e-4b33-ac6a-40a941219854
summer("data/default/CC2_HoldAndKill_n5_t100.csv", t_cc2, sigdigits=4)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
JLD = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
JLD = "~0.13.1"
PlutoUI = "~0.7.38"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Blosc]]
deps = ["Blosc_jll"]
git-tree-sha1 = "310b77648d38c223d947ff3f50f511d08690b8d5"
uuid = "a74b3585-a348-5f62-a45c-50e91977d574"
version = "0.7.3"

[[deps.Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "91d6baa911283650df649d0aea7c28639273ae7b"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.1+0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "63d1e802de0c4882c00aee5cb16f9dd4d6d7c59c"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.1"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.H5Zblosc]]
deps = ["Blosc", "HDF5"]
git-tree-sha1 = "26b22c9039b039e29ec4f4f989946de722e87ab9"
uuid = "c8ec2601-a99c-407f-b158-e79c03c2f5f7"
version = "0.1.0"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires"]
git-tree-sha1 = "e6b1bd8339b2af5a4c2e3103f9dda65de355127e"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.9"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "bab67c0d1c4662d2c4be8c6007751b0b6111de5c"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.1+0"

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

[[deps.JLD]]
deps = ["Compat", "FileIO", "H5Zblosc", "HDF5", "Printf"]
git-tree-sha1 = "958ea85ff8d48a26bf2a3ef4a06ba246fd1b6c7c"
uuid = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
version = "0.13.1"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

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

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═41a85aac-d357-11ec-3d01-71de103a66ec
# ╠═981a2d06-32a8-4fe8-8b10-c9525f72ddf3
# ╠═c4ab5843-067c-44a6-aec9-7f0290f4853a
# ╠═4e207296-a683-4f9c-863b-3205fdda2fbd
# ╠═9665c251-7374-436f-b1fd-c8b5cdef8270
# ╠═9a5471c8-1d00-48ff-a74f-67448835a50c
# ╠═37332697-7b8a-4f6c-a56b-6c29c9ea273c
# ╠═3fbd754c-7ce1-4d2c-83c0-e54fc451ed0e
# ╠═b9e019af-d155-4f1c-9de3-712f1fd02011
# ╠═4fd31c79-c8fe-4d3e-b5a8-e111ae748277
# ╠═d5d11eb5-8f05-4190-aacb-dca923fa894c
# ╠═13ec75c8-874e-4f08-8b8a-416d9f46613c
# ╠═203ebc89-e958-441f-8788-75d1bc945a3f
# ╠═067d97f9-2599-4b64-8791-2d5d9fc49a07
# ╠═db475d07-ed0e-478c-8fb3-5975719266ab
# ╠═f4c0af9a-57bb-4cdc-9f8d-ba6113e59f72
# ╠═a9578595-f100-44dc-89b0-a8926c979acb
# ╠═09ef9182-5350-4347-82e8-fc4cf13c3976
# ╠═596ae945-2d62-4b2a-8349-6c5f416d1dda
# ╠═f1b3a897-70d5-4c68-98f7-6727e72fa442
# ╠═d8029862-0259-4982-94e8-496b4239063e
# ╠═ceb12397-a119-485a-a7e3-fdabfa4ae2f3
# ╠═379c5629-e8df-4188-9716-b0f0db996119
# ╠═692cea24-bf1a-4275-b388-8f2fc7ea610e
# ╠═aba4fefb-f8db-4d37-a287-aa0e6352b4ff
# ╠═edfe4eb1-a007-47d3-a195-704479587365
# ╠═93d654ad-357e-4b33-ac6a-40a941219854
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
