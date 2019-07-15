# combos

This contains the code and processed data for the Chevrier-lab combos project. 

It's components are:
* intermediate data: the processed data used by code in the analysis section

* software: 
  * preprocessing: the scripts for taking raw fastq files into read count tables
    * _pipeline.sh_ - the preprocessing and alignment portion of the ATACseq pipeline. Only handles one set of fastq files at a time and has to be run manually.
      * Data structure: The pipeline is assuming the existence of a project folder where the data is and will be stored. It has a _fastq/_ folder with the fastq files, which are named [stim]_[replicate#]_R[1/2].fastq, and it has an _out/_ folder. _pipe_line.sh_ will populate the _out/_ folder with the results of preprocessing, with the results from each pair of fastq files (each processed by their own call to _pipeline.sh_) being placed in their own folder. The downstream analysis files assume the existence of such folders in _out/_
    * _peak_merge.sh_ - a bash script using bedtools and samtools to make a unified peakset
    * _merge_and_count.r_ - contains an r function that calls peak_merge.r and then makes a table of read counts using the unified peaks
   
  * Analysis: the scripts for running our final analysis, as well as functions they depend on
    * The final folder: each script in this folder runs a portion of our final analysis and produces some number of plots for the paper. Some of the scripts are dependent on having others running for them (and will source the script it depends on when you run it)
      * _final_PI.r_ - analysis of cell proliferation data
        * Figures: 1d, 2a, 2b, 2c, 2d, 2e
        * Supplemental Figures: S4a, S4b, S4c, S5a, S5b, S5c, S6a, S6b, S7, S8a, S8b, S9a, S9b
      * _final_ligand_dilution.r_ - plots of ligand dilution curves
        * Supplemental Figures: S2b
      * _final_SP7.r_ - analysis of mouse RNA-seq data
        * Figures: 3b, 3c, 3d
        * Supplemental Figures: S10, S11
      * _final_human.r_ - analysis of human RNA-seq data
        * Figures: 3e, 3f
        * Supplmental Figures: S12a, S12b, S12c, S12d
      * _final_atac.r_ - analysis of mouse ATAC-seq data
        * Dependencies: _final_SP7.r_ 
        * Supplemental Figures: S13g, S13h
      * _final_secretome.r_ - analysis of mouse secretome mass-spectrometry data and cytokine concentration assays
        * Dependencies: _final_atac.r_
        * Supplemental Figures: S13a, S13b, S13c, S14
      * _final_cyto_restim.R_ - analysis of restimulation assays
        * Figures: 5c, 5d
        * Supplemental Figures: S16a, S17a, S17b
   
    * diff_peaks.r: functions for the differential expression analysis for ATAC-seq data
    * rna_de.r: functions for the differential expression analysis for RNA-seq data
    * combo_analysis.r: Contains a variety of functions usefull for doing analysis of combinations of treatments
      * `compare_combos_contrasts()` runs the bulk of the ATAC-seq data processing
      * `make_comparison_binaries()` calculates which genes are new
      * various other functions for smaller tasks
    * index_functions.r: Contains various functions for Isserlis-formula related analyses
  

