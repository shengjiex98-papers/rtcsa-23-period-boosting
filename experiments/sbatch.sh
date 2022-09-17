#!/bin/bash
#SBATCH --job-name=synthesize-julia
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --mem=4g
#SBATCH --ntasks=8
#SBATCH --time=48:00:00

module add julia
julia --project=.. --threads 8 -e 'import Pkg; Pkg.instantiate(); include("job_dispatch5.jl")'
# julia --project=.. --threads 8 -e 'import Pkg; Pkg.instantiate(); include("experiment.jl"); create_job("CSS", 1, 15, 100, [(0.040, 0.040)], dir="data/common_period_names", clr=false)'
# julia --project=.. --threads 1 -e 'import Pkg; Pkg.instantiate(); include("experiment.jl"); create_job("CC2", 1, 15, 100, [(0.028, 0.028)], dir="data/common_period_names", clr=true, one=(5,6))'
### julia --project=.. --threads 8 experiment.jl "AP" "HK"

### sbatch -p general -N 1 --mem 4096 -n 1 -c 16 -t 48:00:00 --mail-type=end --mail-user=sxunique@cs.unc.edu --wrap='julia --threads 16 experiment.jl "ES" "HK"'