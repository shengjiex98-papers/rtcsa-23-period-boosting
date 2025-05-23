include("emsoft.jl")

function synthesize_traj(bounds, model, strat, t; dims=axes(bounds, 1), hitpattern=repeat([2], t-1))
	# A bug was with `automaton = strat(model[1], model[2], 0, 1)` i.e., it was forbidding any misses to occur.
	automaton = strat(model[1], model[2])
	augbounds = Augment(bounds, automaton)
	# @info bounds, augbounds
	
	get_traj(automaton, augbounds, t; dims=dims, hitpattern=hitpattern)
end

function synthesize(safety_margin, bounds, model, strat, n, max_window_size, t; dims=axes(model[1].A, 1))
	@info "Strat is:" strat
	
	deviations = fill(NaN, (max_window_size-1, max_window_size-1))
	# Windows size starting from 2

	min_hits = 1
	for window_size in 2:max_window_size
		@info "Working on window_size: $(window_size)..."
		while min_hits < window_size
			automaton = strat(model[1], model[2], window_size-min_hits, window_size)
			augbounds = Augment(bounds, automaton)

			iterations = div(t+n-1, n)
			v, i = BoundedTreeIter(automaton, augbounds, n, iterations, safety_margin, dims=dims)
			if v < safety_margin
				# @info "  argmax(d) = $(i)"
				deviations[window_size-1, min_hits] = v
				break
			end

			@info "  maximum(d) = $(v)" window_size min_hits
			min_hits += 1
		end
		@info "window_size: $(window_size) done."
	end
	deviations
end

function synthesize_full(safety_margin, bounds, model, name, strat, n, max_window_size, H; dims=axes(model[1].A, 1), dir="data/default", clr=false)
	@info "Synthesizing for" name strat n H

	deviations = fill(NaN, (max_window_size-1, max_window_size-1))
    indices    = fill(-1,  (max_window_size-1, max_window_size-1))
	time_taken = fill(NaN, (max_window_size-1, max_window_size-1))

	# Windows size starting from 2
	constraints = vcat([[(min_hits, window_size) for min_hits in 1:window_size-1] for window_size in 2:max_window_size]...)

	Threads.@threads for constraint in constraints
		@info "Working on constraint: $(constraint)..."
		min_hits, window_size = constraint

		fullpath = "$(rstrip(dir, '/'))/$(name)_$(bounds[1])_$(min_hits)_$(window_size)_n$(n)_t$(H).csv"
		if !clr && isfile(fullpath)
			# @info "Full path is" fullpath
			v, i, time_elapsed = readdlm(fullpath, ',', Float64)
			@info "Constraint: $(constraint)... loaded from file." (v, i)
		else
			start = time()
			automaton = strat(model[1], model[2], window_size-min_hits, window_size)
			augbounds = Augment(bounds, automaton)
			iterations = div(H+n-1, n)
			v, i = BoundedTreeIter(automaton, augbounds, n, iterations, safety_margin, dims=dims)
			v = round(v, sigdigits=4)
			time_elapsed = round(time() - start, sigdigits=2)
			# @save "data/parallel/$(name)_$(strat)_$(min_hits)_$(window_size).jld" v i time_elapsed
			open(fullpath, "w") do file
				writedlm(file, [v; i; time_elapsed], ',')
			end
			@info "Constraint: $(constraint)... done in $(round(time_elapsed, digits=2)) seconds." (v, i)
		end
		
		deviations[window_size-1, min_hits] = v
		indices[window_size-1, min_hits]    = i
		time_taken[window_size-1, min_hits] = time_elapsed
	end
	
	# @save "data/$(name)_$(strat).jld" deviations indices time_taken
	
	fullpath = "$(rstrip(dir, '/'))/$(name)_$(bounds[1])_n$(n)_t$(H).csv"
	open(fullpath, "w") do file
		writedlm(file, [deviations; indices; time_taken], ',')
	end
	@info "Synthesizing finished for" name strat (deviations, indices, time_taken)
	(deviations, indices, time_taken)
end

function synthesize_one(safety_margin, bounds, model, name, strat, n, min_hits, window_size, t; dims=axes(model[1].A, 1), dir="data/default", clr=false)
	@info "Synthesizing for" name strat n t x_0=bounds[1][1]
	@info "System matrices and K" A=model[1].A B=model[1].B K=model[2]
	@info "Poles are" poles=round.(abs.(eigvals(model[1].A - model[1].B*model[2])), digits=2)'
	@info "Working on constraint: $((min_hits, window_size))..."

	fullpath = "$(rstrip(dir, '/'))/$(name)_$(bounds[1])_$(min_hits)_$(window_size)_n$(n)_t$(t).csv"
	if !clr && isfile(fullpath)
		# @info "Full path is" fullpath
		v, i, time_elapsed = readdlm(fullpath, ',', Float64)
		@info "Constraint: $((min_hits, window_size))... loaded from file." (v, i)
	else
		start = time()
		automaton = strat(model[1], model[2], window_size-min_hits, window_size)
		augbounds = Augment(bounds, automaton)
		iterations = div(t+n-1, n)
		v, i = BoundedTreeIter(automaton, augbounds, n, iterations, safety_margin, dims=dims)
		v = round(v, sigdigits=4)
		time_elapsed = round(time() - start, sigdigits=2)
		# @save "data/parallel/$(name)_$(strat)_$(min_hits)_$(window_size).jld" v i time_elapsed
		open(fullpath, "w") do file
			writedlm(file, [v; i; time_elapsed], ',')
		end
		@info "Constraint: $((min_hits, window_size))... done in $(round(time_elapsed, digits=2)) seconds." (v, i)
	end
	
	@info "Synthesizing finished for" name strat (v, i, time_elapsed)
	(v, i, time_elapsed)
end
