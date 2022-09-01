include("experiment.jl")

periods = [0.020, 0.023, 0.025, 0.027, 0.018, 0.015]
path = "data/closeperiod"

# Recomputed gain values
# for sys_name in sys_names, period in periods
#     create_job(sys_name, 10, (period, period), dir=path, clr=true)
# end

# Old gain values
# for sys_name in sys_names, period in periods
#     create_job(sys_name, 10, (period, 0.020), dir=path, clr=true)
# end

# Same new period, old gain values
# for sys_name in sys_names, period in periods
#     create_job(sys_name, 10, (0.020, period), dir=path, clr=true)
# end

# F1 0.02, 0.02
create_job("F1", 0.1, (0.020, 0.020), dir=path, clr=true)