module Experiment

using ControlTimingSafety

"""
    synthesize_constraints(sysd, K, z_0, d_max, maxwindow, n, t)
Find all `MeetAny` weakly hard constraints with window size at most `maxwindow` that 
guarantees the deviation upper bound is at most `d_max`. The system is specified by 
[`Automaton`](@ref) `a` and initial state is `z_0`. `n` and `t` are as in 
[`bounded_runs_iter`](@ref).
"""
function synthesize_constraints(sysd::AbstractStateSpace{<:Discrete},
    K::AbstractMatrix{Number}, z_0::AbstractVecOrMat, d_max::Float64,
    maxwindow::Integer, n::Integer, t::Integer, nominal_traj::AbstractMatrix{Number})

    safe_constraints = MeetAny[]

    # Do not need to go through all O(maxwindow^2) constraints,
    # see paper for optimization argument
    meet = 1
    for window in 2:maxwindow
        while meet < window
            constraint = MeetAny(meet, window)
            a = hold_kill(sysd, K, constraint)
            # Check if the deviation bound is within the safety margin
            reachable = bounded_runs_iter(a, z_0, n, t)
            m = maximum(deviation(a, z_0, reachable))
            if m <= d_max
                # All constraints with (m, window) where m >= meet are valid
                for i in meet:window-1
                    push!(safe_constraints, MeetAny(i, window))
                end
                break
            end
            meet += 1
        end
    end

    safe_constraints
end

end
