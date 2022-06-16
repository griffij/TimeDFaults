#PBS -P w84
#PBS -q express
#PBS -l ncpus=1
#PBS -l storage=gdata/w84
#PBS -l walltime=00:10:00
#PBS -l mem=4Gb
#PBS -l wd

module load R/4.1.0
cd ../models
R CMD BATCH fit_invgauss_sliprate_ksamples_std_params.R
