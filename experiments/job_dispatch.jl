include("experiment.jl")

periods = [0.005, 0.010, 0.020, 0.040, 0.080, 0.160, 0.320]
path = "data/noewb"

# Recomputed gain values
for sys_name in sys_names, discretization in periods
    create_job(sys_name, 10, (discretization, discretization), dir=path)
end

# Old gain values
for sys_name in sys_names, period in periods
    create_job(sys_name, 10, (discretization, 0.020), dir=path)
end
