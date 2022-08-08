using JLD
include("automata.jl")

sys = let
	A = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
	B = [0.232; 0.0203; 0];
	C = [0 0 1];
	D = [0];
	ss(A, B, C, D)
end
ctrl_delay = 0.1
sysd = c2d(sys, ctrl_delay)
sysd_delay = c2da(sys, ctrl_delay, ctrl_delay)
K = let
	state_cost = 50
	control_cost = 1
	Q = [0 0 0 0; 0 0 0 0; 0 0 state_cost 0; 0 0 0 0]
	R = [control_cost 0; 0 0][1:1,1:1]  # Hack to make a 1x1 matrix.  [1] makes a vector, which has no isposdef method
	lqr(sysd_delay, Q, R)
end
dims = [3]
@info sys
@info sysd
@info sysd_delay
@info K

@save "aircraft_pitch.jld" sys sysd sysd_delay K dims
