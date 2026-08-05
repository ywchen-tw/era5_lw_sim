#!/bin/bash
# CURC SLURM template for the PREFIRE BT simulation + Jacobians (stage 8).
# reptran fine MUST run here, not locally — parallel fine-grid uvspec
# workers exhaust memory on a laptop. Adjust account/partition and REPO.
#
# collocate/prep run locally under the era5 env beforehand; sync the
# derived/YYYY/MM/prefire_bt/ folder (manifest + profile files) to CURC.
#
# Usage:  sbatch slurm/curc_prefire_bt.sh 2025 1 1
#SBATCH --job-name=prefire_bt
#SBATCH --partition=amilan
#SBATCH --account=REPLACE_WITH_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=prefire_bt_%j.out

module purge
module load anaconda intel/2022.1.2 hdf5/1.10.1 zlib/1.2.11 netcdf/4.8.1 swig/4.1.1 gsl/2.7

# er3t env on CURC (the ARCSIX scripts use "er3t"); needs libRadtran + er3t
conda activate er3t

REPO=/projects/$USER/era5_analysis   # <-- adjust
cd "$REPO"

YEAR=${1:?year}; MONTH=${2:?month}; SAT=${3:-1}

python src/prefire_bt.py run --year "$YEAR" --month "$MONTH" --sat "$SAT" \
    --mol-abs-param fine --workers "$SLURM_NTASKS" --overwrite
python src/prefire_bt.py jacobian --year "$YEAR" --month "$MONTH" --sat "$SAT" \
    --simulator lrt --mol-abs-param fine --workers "$SLURM_NTASKS"
