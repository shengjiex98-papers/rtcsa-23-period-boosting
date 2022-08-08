using JLD
include("automata.jl")

sys = let
	R = 0.025
	ωel = 2000*pi
	Ld = 0.0001
	Lq = 0.00012
	A = [-R/Ld  Lq * ωel / Ld;  -Ld * ωel / Lq   -R / Lq]
	B = [1 / Ld   0;  0   1/Lq]
	C = [0 0]
	D = [0 0]
	ss(A, B, C, D)
end
ctrl_delay = 0.00001
sysd = c2d(sys, ctrl_delay)
sysd_delay = c2da(sys, ctrl_delay, ctrl_delay)
K = let
	Q = Array(I(2))
	R = Array(I(2))
	lqr(sysd, Q, R)
end
dims = [1, 2]

@info sys
@info sysd
@info K

@save "electric_steering.jld" sys sysd sysd_delay K dims
