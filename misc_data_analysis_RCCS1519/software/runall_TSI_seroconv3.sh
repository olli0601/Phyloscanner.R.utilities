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
        echo 'qsub -v STEP="xxx" runall_TSI_seroconv3.sh'
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
out_dir_base="$DEEPANALYSES/seroconverters3_alignXX"
out_dir_rel="$out_dir_base/TODO" # TODO
controller="$software_path/runall_TSI_seroconv3.sh" #current script location
CLUSIZE='50'
DATE='2022-07-20'

echo Check that DATE, CLUSIZE and out_dir_rel are well defined.

# helper function to move log files within a directory
# TODO: I think it would be better to have a separate file with functions and source it
move_log_files() {

        DIR=$@
        for dir in $DIR; do

                # skip directories which do not exist
                [  -e "$dir" ] || continue 1 

                for file in $dir/*sh; do
                
                        # skip files that do not exits
                        [ -e "$file" ] || continue 1

                        # create a directory to store all outputs of a single sh script
                        logdir=${file//.sh/}
                        echo "created $DIR/$newdir directory"
                        find ${file}.o* > /dev/null 2>&1 && mkdir -p $logdir
                done
        done
}



cwd=$(pwd)
echo $cwd
module load anaconda3/personal
source activate phylo_alignments

case $STEP in

        # In this analysis we avoid the first step of computing similarities, as there exist already
        net)
        echo "---- initialise analysis ----"
        Rscript $software_path/TSI_initialise_sero2analysis.R \
        --out_dir_base $out_dir_base 
        ;;

        # modified this step adding the reference flag 
        ali)
        echo "---- compute alignments ----"
        Rscript $software_path/make_deep_sequence_alignments.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --windows_start 550 \
        --windows_end 9500 \
        --sliding_width 25 \
        --n_control 0 \
        --cluster_size $CLUSIZE \
        --reference ConsensusGenomes.fasta \
        --mafft "--globalpair --maxiterate 1000" \
        --rm_vloops FALSE \
        --tsi_analysis FALSE
        ;;
        
        # The 2 here should be run without changes the first time
        btr)
        move_log_files $out_dir_base/*_work/
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


        # atm modified from here...
        atr)
        conda activate phylostan
        Rscript $software_path/analyse_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --normRefFileName "$DEEPDATA/normalisation_ByPosition.csv" \
        --outgroupName "REF_CON_H" \
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
        --env_name "phylostan"
        ;;
        
        # ... to here

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

