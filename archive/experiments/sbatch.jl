systems = [
    # Dict([("name", "RC"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    Dict([("name", "F1"), ("x0", 1), ("n", 10), ("p", 0.020), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    # Dict([("name", "DCM"), ("x0", 100), ("n", 10), ("p", 0.023), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    # Dict([("name", "CSS"), ("x0", 100), ("n", 15), ("p", 0.027), ("ctrl", "delay_lqr"), ("ctrl_args", ())]),
    # Dict([("name", "CC2"), ("x0", 1), ("n", 20), ("p", 0.028), ("ctrl", "pole_place"), ("ctrl_args", (0.85))])
]
date_periods = [0.015, 0.028, 0.040]

function setup()
    cmd = `module add julia`
    run(cmd)
end

function build_exec(sys, path, m, k, p, cp)
    # return "'import Pkg; Pkg.instantiate(); include(\"experiment.jl\"); create_job(\"$(sys["name"])\", $(sys["x0"]), $(sys["n"]), 100, ($p, $cp); dir=\"$path\", clr=true, ctrl=\"$(sys["ctrl"])\", one=($m, $k))'"
    return "'include(\"experiment.jl\"); create_job(\"$(sys["name"])\", $(sys["x0"]), $(sys["n"]), 100, ($p, $cp); " *
            "dir=\"$path\", clr=true, ctrl=\"$(sys["ctrl"])\", ctrl_args=$(sys["ctrl_args"]), one=($m, $k))'"
end

function releasejob(j_exec)
    wrap = "julia --project=.. -e $j_exec"
    # j_program = "julia"
    # j_options = ["--project=..", "-e $j_exec"]
    # wrap = "$j_program $j_options"

    program = `sbatch -p general -N 1 -n 1 --mem=2g -t 00-48:00:00 --wrap=$wrap`
    # program = "sbatch"
    # options = ["-p general", "-N 1", "-n 1", "--mem=2g", "-t 00-48:00:00", "--wrap=\"$wrap\""]
    print(program, '\n')
    run(program)
end

# setup()
for sys in systems, p in date_periods, cp in [sys["p"], p], k in 2:6, m in 1:k-1
    releasejob(build_exec(sys, "data/ccs_lqr/", m, k, p, cp))
    # break
    # sleep(1)
end