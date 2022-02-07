module purge

export TMPDIR=${PBS_JOBFS:-/tmp}

function useconda() {
	                eval "$(/g/data/xf3/ht5438/miniconda3/bin/conda shell.zsh hook)"
			        }


useconda

source /g/data/xf3/ht5438/miniconda3/etc/profile.d/conda.sh
conda activate snakemake