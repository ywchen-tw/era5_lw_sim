#!/bin/bash
# CURC SLURM template for the LW flux simulation (stage 7).
# Modeled on ~/programming/arcsix_sfc/lrt_sim/curc_shell_*.sh — adjust the
# account/partition and the REPO path to your CURC layout before submitting.
#
# Usage:  sbatch slurm/curc_lrt_sim.sh 2020 1 1 12
#SBATCH --job-name=lrt_sim
#SBATCH --partition=amilan
#SBATCH --account=REPLACE_WITH_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=04:00:00
#SBATCH --output=lrt_sim_%j.out

module purge
module load anaconda intel/2022.1.2 hdf5/1.10.1 zlib/1.2.11 netcdf/4.8.1 swig/4.1.1 gsl/2.7

# er3t env on CURC (the ARCSIX scripts use "er3t"); needs libRadtran + er3t
conda activate er3t

REPO=/projects/$USER/era5_analysis   # <-- adjust
cd "$REPO"

YEAR=${1:?year}; MONTH=${2:?month}; DAY=${3:?day}; HOUR=${4:?hour}

python src/lrt_sim.py prep --year "$YEAR" --month "$MONTH" --day "$DAY" --hour "$HOUR"
python src/lrt_sim.py run  --year "$YEAR" --month "$MONTH" --day "$DAY" --hour "$HOUR" \
    --streams 8 --mol-abs-param fine --workers "$SLURM_NTASKS"
python src/lrt_sim.py compare --year "$YEAR" --month "$MONTH" --day "$DAY" --hour "$HOUR"
