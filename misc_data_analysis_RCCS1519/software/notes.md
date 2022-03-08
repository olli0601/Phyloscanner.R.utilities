# Phyloscanner Pipeline for TSI Analysis

Here I write my comments and notes on XX code. Objectives:

1. Make sure this can be run in a pipeline
2. Make sure I have the environments to run this on the HPC
3. Understand the specific components: what each part is doing + what data are required + where are the results stored
4. Try to think about how I can generalise this to other settings.



**QUESTIONS:**

- As some steps in the analysis are extremely expensive, we want to be careful in many aspects of our analysis: Who are the eligible individuals for analysis?

  - ATM people with consensus sequences (eg in here)
  - But, to run Tanya, we need at least one basefreq file for each participant.

  Can we check that we have those? If not, ask Oli/Tanya on what to do.
  
- Also, as a feasibility check, why don t we just train something similar to Tanya's algorithm by ourselves? They say linear regression is not too bad, so why not just extract predictors we are interested in and run lin reg? (this may not be easy, but still should not take weeks as TSI analysis)





- OLI: do we want excised positions or not?
- XX: can we have an option not to excise them?
- do not really like how sometimes out.dir and args$out.dir point to different locations. What to do about this



**TODO:**

- need to pass the --infile from TSIcalculate_similarity_over_genome.R to calculate_similarity (shouldn t be hard)
- ldpaths issue (see [here](https://github.com/conda-forge/r-base-feedstock/issues/67)): for the moment I use a while loop that re-runs the Rscript until the output is produced. However, the conda env was called at the beginning of the script. Should I activate the environment justbefore running the script?





FILL THIS TABLE IN AS YOU WORK THROUGH MRC EXAMPLE

| SCRIPT                                | INPUT                                                        | OUTPUT                                                       | env                   | step |
| ------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | --------------------- | ---- |
| `TSIcalculate_similarity_over_genome` | [consensus seqs](/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_MRC_alignment.fasta) | [bash scripts](script_calculate_similarity_job1.qsub)        | [1](phylo_alignments) | sim  |
| ↘ `calculate_similarity`<br />        |                                                              | [similarity scores](out_dir_base/potential_network/similarityX.rds) | [1](phylo_alignments) |      |
| `TSIcreate_network`                   | similarity scores +<br />[db_sharing_extract](PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv ---- PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv' )<br />[phscinput_samples]('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phsci    nput_samples.rds')<br /> | [cluster assignments](clusters.rds) <br />[phscinput](phscinput_runs_clusize_100_ncontrol_0.rds) | [1](phylo_alignments) | net  |
| `make_deep_sequence_alignments`       | above                                                        | bash scripts 4<br />[groupings](out_dir_base/potential_network/phscinput_runs_clusize_X_ncontrol_Y.rds.) | [1](phylo_alignments) |      |
| BUILD TREE                            |                                                              |                                                              |                       |      |
| CHECK TREES                           |                                                              |                                                              |                       |      |
| `analyse_trees.R`                     |                                                              |                                                              |                       |      |
|                                       |                                                              |                                                              |                       |      |
|                                       |                                                              |                                                              |                       |      |



**Thinking about .sh controller**

Aim is to direct the full analysis with a unique `.sh` script. The problem is that we cannot run one script after the other, as each script produces some child jobs that may take loads of time to finish.
The idea, is then to have the last of this completed jobs to 'inform' the controller that every analysis is done. I have to refine the idea, but below are useful facts:

- `qsub -v VAR1="value1",VAR2="value2" job.sh` allows to pass variables inside the scripts.
  This could be useful to inform the controller on what point we are in the analysis. Found [here](https://stackoverflow.com/questions/18925068/how-to-pass-parameters-from-qsub-to-bash-script)
- at the moment working on it at `~.scripts/testing` on HPC.
- can call this something like: `runall_TSI.sh`



**MEETING 2022 02 28**

- asked few questions to XX . No surprising responses. Come back to see if she fixed the `html`
- 



**INTERESTING**

- 92% of inds in the 200422_PANGEA2_RCCS_alignment.fasta files have a corresponding base frequency file in the Fraser directories. (680 don't and) (MOST ARE PG19, 40 are PG15).
- 98% of inds in the 200422_PANGEA2_MRC_alignment.fasta files have a corresponding base frequency file in the Fraser directories (36 do not, all are PG15).



