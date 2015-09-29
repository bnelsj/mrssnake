module purge

module load modules modules-init modules-gs/prod modules-eichler/prod 

module use /net/eichler/vol7/home/psudmant/local_modules
export MOD_PS_SW=/net/eichler/vol7/home/psudmant/local_modules_sw

module load zlib/1.2.6 mrsfast/2.5.0.4_READ_PATCHED
module load samtools/1.2
module load anaconda/2.3.0
