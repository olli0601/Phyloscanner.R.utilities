#!/bin/sh
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -lwalltime=05:00:00
#PBS -j oe 

# The key driver of this analysis is the STEP parameter
# which should be passed through the qsub command.
# qsub -v STEP="net" runall_TSI_seroconv2.sh
# If unset default to "sim"
if [ -z "$STEP" ]
then
        echo "Intended use:\n"
        echo 'qsub -v STEP="xxx" runall_TSI.sh'
        exit 1
fi
echo "running '${STEP:=sim}' analysis"

# This includes all code necessary to run PHSC pipeline to produce TSI estimates
DEEPDATA="/rds/general/project/ratmann_pangea_deepsequencedata/live"
DEEPANALYSES="/rds/general/project/ratmann_deepseq_analyses/live"
HOME="/rds/general/user/ab1820/home"

software_path="$HOME/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software"
phyloscanner_path="$HOME/git/phyloscanner"
hivtsipath="$HOME/git/HIV-phyloTSI"
out_dir_base="$DEEPANALYSES/PANGEA2_RCCS1519_UVRI_TSI"
out_dir_rel="$out_dir_base" # TODO
controller="$software_path/runall_TSI.sh" 
CLUSIZE='50'
DATE='2022-05-19'


cwd=$(pwd)
echo $cwd
module load anaconda3/personal
source activate phylo_alignments

case $STEP in

        net)
        echo "---- initialise analysis ----"
        Rscript TSI_initialise.R \
        --out_dir_base $out_dir_base \
        --cluster_size $CLUSIZE \
        --controller $controller
        ;;

        ali)
        echo "---- compute alignments ----"
        Rscript $software_path/make_deep_sequence_alignments.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --sliding_width 10 \
        --n_control 0 \
        --cluster_size $CLUSIZE \
        --tsi_analysis TRUE
        ;;
        
        btr)
        echo "----- build trees ----"
        Rscript $software_path/make_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE 
        ;;

        ctr)
        echo "----- check trees ----"
        Rscript $software_path/check_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE
        ;;

        atr)
        conda activate phylostan
        Rscript $software_path/analyse_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --blacklistReport \
        --outputNexusTree \
        --skipSummaryGraph \
        --normRefFileName "$phyloscanner_path/InfoAndInputs/HIV_DistanceNormalisationOverGenome.csv" \
        --outgroupName "A1.UGANDA.2007.p191845.JX236671"  \
        --ratioBlacklistThreshold 0.005 \
        --date $DATE \
        --env_name "phylostan"
        ;;

        tsi)
        echo "----- Run HIV-TSI -----"
        Rscript $software_path/TSI_run_predictions.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --input_samples "$out_dir_base/220419_phscinput_samples.rds" \
        --TSI_dir $hivtsipath \
        --date $DATE \
        --env_name 'hivphylotsi'
        ;;

        dti)
        echo "----- get dates of infection -----"
        Rscript $software_path/TSI_estimate_dates.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --date $DATE \
        --input_samples "$out_dir_base/220419_phscinput_samples.rds" 
        ;;

        *)
        echo "no R script run. STEP does not match any task.\n" 
        ;;
esac
