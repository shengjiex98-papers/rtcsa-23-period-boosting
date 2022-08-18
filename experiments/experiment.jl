include("iccad.jl")

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

model_map = Dict(
    "DCM" => (c2d(sys_dcm, h), delay_lqr(sys_dcm, h)),
    "CSS" => (c2d(sys_css, h), delay_lqr(sys_css, h)),
    "EWB" => (c2d(sys_ewb, h), delay_lqr(sys_ewb, h)),
    "CC1" => (c2d(sys_cc1, h), delay_lqr(sys_cc1, h)),
    "CC2" => (c2d(sys_cc2, h), delay_lqr(sys_cc2, h)),
    "RC" => (c2d(sys_rc, h), delay_lqr(sys_rc, h)),
    "F1" => (sysd_f1, [0.293511 0.440267])
)
model_names = sort([keys(model_map)...])

function create_job(model_name, s0, strat_names...; dir="data/default", clr=false)
    @info "Threads: " Threads.nthreads()
    model = model_map[model_name]

    n = 5
    t = 100
    max_window_size = 6
    safety_margin = 1000

    q = size(model[1].A, 1)
    bounds = repeat([s0 s0], q)

    for strat_name in strat_names
        strat = strat_map[strat_name]
        synthesize_full(safety_margin, bounds, model, model_name, strat, n, max_window_size, t; dims=[2], dir=dir, clr=clr)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Running as main"
    create_job(ARGS[1], ARGS[2:end]...)
else
    @info "Included as module"
end