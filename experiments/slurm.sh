#!/bin/bash
#SBATCH --job-name=synthesize-julia
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --mem=4g
#SBATCH --ntasks=8
#SBATCH --time=48:00:00

module add julia
julia --project=.. --threads 8 -e 'import Pkg; Pkg.instantiate(); include("job_dispatch.jl")'
### julia --project=.. --threads 8 experiment.jl "AP" "HK"

### sbatch -p general -N 1 --mem 4096 -n 1 -c 16 -t 48:00:00 --mail-type=end --mail-user=sxunique@cs.unc.edu --wrap='julia --threads 16 experiment.jl "ES" "HK"'