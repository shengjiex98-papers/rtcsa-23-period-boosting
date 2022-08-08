using JLD
include("automata.jl")

sys = let
	A = [1 0.13; 0 1];
	B = [0.02559055118110236; 0.39370078740157477];
	C = [0 0];
	D = [0];
	ss(A, B, C, D)
end
ctrl_delay = 0.02
sysd = let
	A = [1 0.13; 0 1];
	B = [0.02559055118110236; 0.39370078740157477];
	C = [0 0];
	D = [0];
	ss(A, B, C, D, ctrl_delay)
end
sysd_delay = c2da(sys, ctrl_delay, ctrl_delay)
K = [0.293511 0.440267]
dims = [1, 2]
@info sys
@info sysd
@info sysd_delay
@info K

@save "f1tenth_car.jld" sys sysd sysd_delay K dims
