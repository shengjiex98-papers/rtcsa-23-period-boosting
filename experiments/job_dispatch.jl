include("experiment.jl")
using Combinatorics

execution_times = Dict(
    "RC" => 0.010,
    "F1" => 0.013,
    "DCM" => 0.012,
    "CSS" => 0.010,
    "CC2" => 0.015
)

# function common_period(execution_times, num_controllers)
function common_period(execution_times, num_controllers=2)
    return Set(sum.(collect(combinations(collect(values(execution_times)), num_controllers))))
end

@info common_period(execution_times)

# periods = [0.020, 0.023, 0.025, 0.027, 0.018, 0.015]
periods = common_period(execution_times)
path = "data/common_period_names"

sys_names = ["RC"]
x_0 = 100

# Recomputed gain values
for sys_name in sys_names, period in periods
    create_job(sys_name, x_0, (period, period), dir=path, clr=false)
    create_job(sys_name, x_0, (period, 0.028), dir=path, clr=false)
    create_job(sys_name, x_0, (0.028, period), dir=path, clr=false)
end