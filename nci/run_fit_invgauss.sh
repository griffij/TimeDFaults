#PBS -P w84
#PBS -q express
#PBS -l ncpus=1
#PBS -l storage=gdata/w84
#PBS -l walltime=02:00:00
#PBS -l mem=8Gb
#PBS -l jobfs=2Gb
#PBS -l wd


module load R/4.1.0
export LD_LIBRARY_PATH=/home/547/jdg547/lib/:$LD_LIBRARY_PATH

cd ../models
R CMD BATCH fit_invgauss_sliprate_ksamples_std_params_hyde.R
