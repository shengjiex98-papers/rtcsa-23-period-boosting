include("experiment.jl")
using Combinatorics

execution_times = Dict(
    "RC" => 0.010,
    "F1" => 0.013,
    "DCM" => 0.012,
    "CSS" => 0.010,
    "CC2" => 0.015
)

systems = [
    # Dict([("name", "RC"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", delay_lqr)]),
    # Dict([("name", "F1"), ("x0", 1), ("n", 16), ("p", 0.020), ("ctrl", delay_lqr)]),
    # Dict([("name", "DCM"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", delay_lqr)]),
    # Dict([("name", "CSS"), ("x0", 100), ("n", 16), ("p", 0.027), ("ctrl", delay_lqr)]),
    Dict([("name", "CC2"), ("x0", 1), ("n", 20), ("p", 0.028), ("ctrl", pole_place)])
]

# function common_period(execution_times, num_controllers)
function common_period(execution_times, num_controllers=2)
    return Set(round.(sum.(collect(combinations(collect(values(execution_times)), num_controllers))), sigdigits=2))
end

@info common_period(execution_times)

# periods = [0.020, 0.023, 0.025, 0.027, 0.018, 0.015]
# periods = common_period(execution_times)
path = "data/finalize"
date_periods = [0.015, 0.028, 0.040]

# Each system, independently
for sys in systems
    create_job(sys["name"], sys["x0"], sys["n"], 100, (sys["p"], sys["p"]); dir=path, clr=true, ctrl=sys["ctrl"])
end

# Each system, using each common period
for sys in systems, period in date_periods
    create_job(sys["name"], sys["x0"], sys["n"], 100, (period, sys["p"]); dir=path, clr=true, ctrl=sys["ctrl"])
    create_job(sys["name"], sys["x0"], sys["n"], 100, (period, period); dir=path, clr=true, ctrl=sys["ctrl"])
end

# sys_names = ["CC2"]
# x_0 = 100

# create_job("F1", 1, (0.028, 0.020), dir=path, clr=false)
# create_job("CC2", 1, 16, 100, [(0.028, 0.028)], dir=path, clr=false, one=[5, 6])
# create_job("CSS", 100, (0.028, 0.02), dir=path, clr=false)

# Recomputed gain values
# for sys_name in sys_names, period in periods
#     create_job(sys_name, x_0, (period, period), dir=path, clr=false)
#     create_job(sys_name, x_0, (period, 0.028), dir=path, clr=false)
#     create_job(sys_name, x_0, (0.028, period), dir=path, clr=false)
# end