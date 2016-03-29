# mrssnake
Simple mrsfast read depth mapping using snakemake

## Quick start
1. Clone the repository, set up the config file

   ```bash
   git clone --recursive https://github.com/bnelsj/mrssnake
   ```
2. Create a tab-delimited manifest file with the appropriate header and a line for each sample

   | sn | source | bam | index |
   | --- | ------ | --- | ----- |
   | sample_name  |  local | /full/path/to/bam | /full/path/to/bam/index |

3. Modify `config.yaml`

   In particular, set the manifest variable to point to your manifest file, 
   make sure the `reference` variable points to the appropriate reference, 
   and that the paths for the `masked_ref` and `contigs` for your reference are correct.
   
   To add a new reference `new_ref` to the config file, change `reference`:
   ```yaml
   reference: new_ref
   ```
   And add lines for `new_ref`:
   ```yaml
   new_ref:
       masked_ref: /path/to/masked_ref
       contigs: /path/to/contigs
   ```
   Note that every `masked_ref` must be TRF and repeat-masked and indexed with mrsfastULTRA.
   
4. Run snakemake

   This example will use 100 cores and the `--cluster-sync` option:
   ```python
   snakemake -T --cluster-sync "qsub -sync y {params.sge_opts}" -w 30 -j 100
   ```


