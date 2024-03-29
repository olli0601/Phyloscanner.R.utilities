#!/bin/sh
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -lwalltime=05:00:00
#PBS -j oe 

# The key driver of this analysis is the STEP parameter
# which should be passed through the qsub command.
# qsub -v STEP="net" runall_TSI_seroconv.sh
# If unset default to "sim"
if [ -z "$STEP" ]
then
        echo "Intended use:\n"
        echo 'qsub -v STEP="xxx" runall_TSI_seroconv.sh'
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
out_dir_base="$DEEPANALYSES/seroconverters"
out_dir_rel="$out_dir_base/2022_03_31_phsc_phscrelationships_sd_42_blacklist_report_T_mr_1_og_REF_CON_H_output_nexus_tree_T_rtt_0005_skip_summary_graph_T_sdt_1/"
controller="$software_path/runall_TSI_seroconv.sh" #current script location
CLUSIZE='6'



cwd=$(pwd)
echo $cwd
module load anaconda3/personal
source activate phylo_alignments

case $STEP in

        # In this analysis we avoid the first step of computing similarities, as there exist already
        net)
        echo "---- get closest individuals ----"
        Rscript $software_path/TSI_find_close_to_sero.R \
        --cluster_size $CLUSIZE \
        --out_dir_base $out_dir_base \
        --sero_tsi "$DEEPDATA/PANGEA2_RCCS/220329_TSI_seroconverters.csv" \
        --path_similarities "$DEEPANALYSES/test_full/potential_network/whole_genome_average_similarity.rds"\
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
        --cluster_size $CLUSIZE
        ;;
        
        btr)
        echo "----- build trees ----"
        Rscript $software_path/make_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date "2022-03-31"
        ;;

        ctr)
        echo "----- check trees ----"
        Rscript $software_path/check_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date "2022-03-31" 
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
        --outgroupName "REF_CON_H"  \
        --ratioBlacklistThreshold 0.005 \
        --date "2022-03-31" \
        --env_name "phylostan"
        ;;

        tsi)
        echo "----- Run HIV-TSI -----"
        Rscript $software_path/TSI_run_predictions.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --input_samples "$out_dir_base/210120_RCCSUVRI_phscinput_samples.rds" \
        --TSI_dir $hivtsipath \
        --date "2022-03-31" \
        --env_name 'hivphylotsi'
        ;;

        dti)
        echo "----- get dates of infection -----"
        Rscript $software_path/TSI_estimate_dates.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --date "2022-03-31" \
        --input_samples "$out_dir_base/210120_RCCSUVRI_phscinput_samples.rds"
        ;;

        pst)
        conda activate phylostan
        echo "----- plot results -----"
        Rscript $software_path/TSI_postprocessing_comparison.R \
        --relationship_dir $out_dir_rel \
        --TSI_dir $hivtsipath \
        --input_samples "$out_dir_base/210120_RCCSUVRI_phscinput_samples.rds"
        ;;

        *)
        echo "no R script run. STEP does not match any task." 
        ;;
esac

