ran pipeline using the command:

conda activate snakemake
module load Singularity

snakemake --profile prf_PanGenotyping/ --singularity-args "--bind /gpfs"
