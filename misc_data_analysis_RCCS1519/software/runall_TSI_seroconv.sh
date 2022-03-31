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
#out_dir_out="$out_dir_base/19037_phsc_output"
#out_dir_rel="$out_dir_base/19037_phsc_phscrelationships_seed_42_blacklist_report_TRUE_distance_threshold_1_min_reads_per_host_1_multinomial_TRUE_outgroup_name_BFR83HXB2_LAI_IIIB_BRUK03455_output_nexus_tree_TRUE_ratio_blacklist_threshold_0005_skip_summary_graph_TRUE/"
#out_dir_logs="$HOME/.scripts/testing"
controller=$(readlink -f ${BASH_SOURCE}) #current script location
CLUSIZE='6L'



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
        --date "2022-03-31"
        ;;

        ctr)
        echo "----- check trees ----"
        Rscript $software_path/check_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --date "2022-03-31"
        ;;
        # Steps here to implement as we go.

        tsi)
        echo "----- Run HIV-TSI -----"
        Rscript $software_path/TSI_run_predictions.R \
        --out_dir_base $out_dir_out \
        --relationship_dir $out_dir_rel \
        --input_samples "$out_dir_base/210120_RCCSUVRI_phscinput_samples.rds" \
        --TSI_dir $hivtsipath \
        --env_name 'hivphylotsi'
        ;;

        dti)
        echo "----- get dates of infection -----"
        Rscript $software_path/TSI_estimate_dates.R \
        --out_dir_base $out_dir_out \
        --relationship_dir $out_dir_rel \
        --input_samples "$out_dir_base/210120_RCCSUVRI_phscinput_samples.rds"
        ;;

        *)
        echo "no R script run. STEP does not match any task." 
        ;;
esac

