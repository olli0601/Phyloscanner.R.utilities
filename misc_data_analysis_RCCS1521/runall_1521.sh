#!/bin/sh
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -lwalltime=05:00:00
#PBS -j oe 

# The key driver of this analysis is the STEP parameter
# which should be passed through the qsub command.
# qsub -v STEP="net" runall_TSI_seroconv2.sh
# If unset default to "none"
if [ -z "$STEP" ]
then
        echo "Intended use:\n"
        echo 'qsub -v STEP="xxx",RES=0 runall_1521.sh'
        echo "OPTIONs:"
        echo "    STEP : one of ali, btr, atr                           [default: none]"
        echo "    RES: determines resources fr pbs jobs (1 to 3)        [default: 1]"
        echo "    N_JOB_RUNNER: number of people that can submit jobs   [default: 1]"
        exit 1
fi

${RES:=1} 
${REDO:=0}
${N_JOB_RUNNER:=1}
echo "running '${STEP:=none}' analysis"

# This includes all code necessary to run PHSC pipeline to produce TSI estimates
DEEPDATA="/rds/general/project/ratmann_pangea_deepsequencedata/live"
DEEPANALYSES="/rds/general/project/ratmann_deepseq_analyses/live"
HOME="/rds/general/user/ab1820/home"

software_path="$HOME/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521"
input_samples="$out_dir_base/240809_RCCSUVRI_phscinput_runs.rds"
phyloscanner_path="$HOME/git/phyloscanner"
hivtsipath="$HOME/git/HIV-phyloTSI"
out_dir_base="$DEEPANALYSES/PANGEA2_RCCS1521"
# out_dir_rel="$out_dir_base/..."
controller="$software_path/$PBS_JOBNAME" #current script location
# For TSI's 25 is fine, for networks 10:
sliding_width=10

CLUSIZE='50'
DATE='2024-09-20'

echo Check that DATE, CLUSIZE and out_dir_rel are well defined.

cwd=$(pwd)
echo $cwd
module load anaconda3/personal
source activate phylo_alignments

case $STEP in

    # The initial steps producing nets were done by Xiaoyue

        # modified this step adding the reference flag 
        # This also has walltime flag but don't know what to do exactly about it...
        ali)
        echo "---- compute alignments ----"

        if [ "$REDO" = "0"]; then
                Rscript $software_path/Rk1521_03_make_deep_sequence_alignments.R \
                --out_dir_base $out_dir_base \
                --pkg_dir $software_path \
                --prog_dir $phyloscanner_path \
                --windows_start 550 \
                --windows_end 9500 \
                --sliding_width $sliding_width \
                --n_control 0 \
                --cluster_size $CLUSIZE \
                --reference ConsensusGenomes.fasta \
                --mafft "--globalpair --maxiterate 1000 " \
                --rm_vloops \
                --controller $controller \
                --walltime_idx $RES \
                --tsi_analysis FALSE \
                --split_jobs_by_n $N_JOB_RUNNER
        else
                Rscript $software_path/Rk1521_03_make_deep_sequence_alignments.R \
                --out_dir_base $out_dir_base \
                --pkg_dir $software_path \
                --prog_dir $phyloscanner_path \
                --windows_start 550 \
                --windows_end 9500 \
                --sliding_width $sliding_width \
                --n_control 0 \
                --cluster_size $CLUSIZE \
                --reference ConsensusGenomes.fasta \
                --mafft " --globalpair --maxiterate 1000 " \
                --rm_vloops \
                --controller $controller \
                --walltime_idx $RES \
                --date $DATE \
                --tsi_analysis FALSE \
                --split_jobs_by_n $N_JOB_RUNNER
        fi
        ;;
        
        # The 2 here should be run without changes the first time
        # I believe there is no reason for having 2 separate scripts...
        btr)
        echo "----- build trees ----"
        Rscript $software_path/Rk1521_04_make_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE \
        --controller $controller \
        --walltime_idx $RES
        ;;

        # We may not need this.
        ctr)
        echo "----- check trees ----"
        Rscript $software_path/Rk1521_05_check_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE \
        --walltime_idx $RES
        ;;

        # DOUBLE CHECK HERE AGAINST ORIGINAL!!!
        atr)
        conda activate phylostan
        Rscript $software_path/Rk1521_06_analyse_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --normRefFileName "$DEEPDATA/normalisation_ByPosition.csv" \
        --outgroupName "REF_B.FR.83.HXB2_LAI_IIIB_BRU_K03455" \
        --ratioBlacklistThreshold 0.01 \
        --distanceThreshold "0.02 0.05"   \
        --maxReadsPerHost 100 \
        --minReadsPerHost 30  \
        --multinomial \
        --noProgressBars  \
        --postHocCountBlacklisting  \
        --relaxedAncestry \
        --zeroLengthAdjustment \
        --date $DATE \
        --controller $controller \
        --env_name "phylostan" \
        --verbose TRUE
        ;;

        *)
        echo "no R script run. STEP does not match any task.\n" 
        ;;
esac

