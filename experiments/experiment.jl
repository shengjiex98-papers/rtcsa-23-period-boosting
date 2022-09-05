include("aspdac.jl")

sys_rc = let
	r_1 = 100000
	r_2 = 500000
	r_3 = 200000
	c_1 = 0.000002
	c_2 = 0.000010
	A = [-1/c_1 * (1/r_1 + 1/r_2)  1/(r_2*c_1)
	     1/(r_2*c_2)               -1/c_2 * (1/r_2 + 1/r_3)]
	#B = [1/(r_1*c_1)  0
	#     0            1/(r_3*c_2)]
	B = [1/(r_1*c_1)
	     1/(r_3*c_2)]
	C = [1  -1]
	D = 0

	ss(A, B, C, D)
end

# This is discretized with period h=0.02s
sysd_f1 = let
	A = [1 0.13; 0 1];
	B = [0.02559055118110236; 0.39370078740157477];
	C = [0 0];
	D = [0];
	ss(A, B, C, D, 0.02)
end

sys_dcm = let
    A = [-10 1; -0.02 -2]
    B = [0; 2]
    C = [1 0]
    D = 0

    ss(A, B, C, D)
end

sys_css = let
    A = [0 1 0 0; -8 -4 8 4; 0 0 0 1; 80 40 -160 -60]
    B = [0; 80; 20; -1120]
    C = [1 0 0 0]
    D = 0

    ss(A, B, C, D)
end

sys_ewb = let
    A = [0 1; 8.3951e3 0]
    B = [0; 4.0451]
    C = [7.9920e3 0]
    D = 0

    ss(A, B, C, D)
end

sys_cc1 = let
    A = -0.05
    B = 0.01
    C = 1
    D = 0

    ss(A, B, C, D)
end

sys_cc2 = let
    A = [0 1 0; 0 0 1; -6.0476 -5.2856 -0.238]
    B = [0; 0; 2.4767]
    C = [1 0 0]
    D = 0

    ss(A, B, C, D)
end

h = 0.02

delay_lqr(sys, h) = let
    Q = I
    R = I
    sysd_delay = c2da(sys, h, h)
    lqr(sysd_delay, Q, R)
end

# sys_map = Dict(
#     "RC" => (c2d(sys_rc, h), delay_lqr(sys_rc, h)),
#     # "F1" => (sysd_f1, [0.293511 0.440267]),
#     "DCM" => (c2d(sys_dcm, h), delay_lqr(sys_dcm, h)),
#     "CSS" => (c2d(sys_css, h), delay_lqr(sys_css, h)),
#     "EWB" => (c2d(sys_ewb, h), delay_lqr(sys_ewb, h)),
#     "CC1" => (c2d(sys_cc1, h), delay_lqr(sys_cc1, h)),
#     "CC2" => (c2d(sys_cc2, h), delay_lqr(sys_cc2, h))
# )
sys_map = Dict(
    "RC" => sys_rc,
    "DCM" => sys_dcm,
    "CSS" => sys_css,
    # "EWB" => sys_ewb,
    "CC1" => sys_cc1,
    "CC2" => sys_cc2
)
sys_names = sort([keys(sys_map)...])

function create_job(sys_name, x0, periods...; dir="data/default", clr=false)
    @info "Threads: " Threads.nthreads()

    if !isdir(dir)
        mkdir(dir)
    end

    system = sys_map[sys_name]

    n = 5
    t = 100
    max_window_size = 6
    safety_margin = 1000

    q = size(system.A, 1)
    bounds = repeat([x0 x0], q)

    for hs in periods
        @info "Period is" hs[1] hs[2]
        subdir = "$(rstrip(dir, '/'))/$(hs[1])s_$(hs[2])s"
        if !isdir(subdir)
            mkdir(subdir)
        end

        strat = HoldAndKill
        # model = (sysd_f1, [0.293511 0.440267])
        model = (c2d(system, hs[1]), delay_lqr(system, hs[2]))
        
        synthesize_full(safety_margin, bounds, model, sys_name, strat, n, max_window_size, t; dims=[2], dir=subdir, clr=clr)
    end

    @info "Experiment finished."
end

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Running as main"
    # create_job(ARGS[1], ARGS[2:end]...)
    # print(ARGS)
else
    @info "Included as module"
end