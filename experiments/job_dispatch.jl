include("experiment.jl")
using Combinatorics

# sys_map = Dict(
#     "RC" => sys_rc,
#     "DCM" => sys_dcm,
#     "CSS" => sys_css,
#     # "EWB" => sys_ewb,
#     "CC1" => sys_cc1,
#     "CC2" => sys_cc2
# )

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
path = "data/common_period"

sys_names = ["F1"]

# Recomputed gain values
for sys_name in sys_names, period in periods
    create_job(sys_name, 1, (period, period), dir=path, clr=true)
end

# Old gain values
for sys_name in sys_names, period in periods
    create_job(sys_name, 1, (period, 0.028), dir=path, clr=true)
end

# Same new period, old gain values
for sys_name in sys_names, period in periods
    create_job(sys_name, 1, (0.028, period), dir=path, clr=true)
end

# # F1 0.02, 0.02
# create_job("F1", 0.1, (0.020, 0.020), dir=path, clr=true)