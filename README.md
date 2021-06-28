# STXBP1

### Requirements:
* [R](https://www.r-project.org/) with packages tidyverse, stringr, dplyr, Hmisc, memoise, reshape2, readr, corrplot, and RColorBrewer.

### Steps to Run:
* Clone the repository, modify the [config](https://github.com/helbig-lab/STXBP1/blob/master/input.yml) file as needed.

* In the [config file](https://github.com/helbig-lab/SCN2A/blob/master/input.yml) determine the main output_dir, this is where your output files would be written to.  Default parameters have been added for most aspects of the [config file](https://github.com/helbig-lab/STXBP1/blob/master/input.yml). Processing methods and algorithms of calculation for similarity analysis can be specified by the user in the terminal after running the below script. 

* Run [R file](https://github.com/helbig-lab/STXBP1/blob/master/master_config.R), specifying the YAML config file using the --input flag.

```
~/Rscript master_config.R --input /path_to/input.yml
```

### Running the tests
There are test files available here: [Files](https://github.com/helbig-lab/STXBP1/tree/master/raw_files). Ensure that these files are linked appropriately in the [config file](https://github.com/helbig-lab/STXBP1/blob/master/input.yml) as such:

```

file_path : raw_files/STXBP1_full_base_v.csv

```

This provides the necessary cohort of patients with annotated HPO terms.

Note that there is an option to manually perform term propagation (see below), though default files have been provided.

## Phenotypic Similarity Analysis
Using the Human Phenotype Ontology (HPO) and a cohort of individuals annotated via HPO terms and VCF files, these scripts find phenotype-genotype correlations. This can aid in gene discovery, treatment, and a better understanding of genes' phenotypic variability. First, HPO terms' similarity scores are found for every individual pair using either the Resnik or Cube algorithm. Next, genes with potentially causitive variants in multiple individuals are extracted and the median similarity scores among each of these patients is calculated. Using permutation analysis of median similarity scores, p-values are assigned to each of these genes. A lower p-value (p < 0.05) potentially indicates a causitive gene.

### Warning
* Running the similarity analyses scripts outside of a cluster or powerful computer may take an excessively long time.  It is recommended that you run the similarity analyses in a cluster or reduce total data in your copy of the example variant file to speed up the process.

* If the ```gene_count_cube_auto.R``` file does not run, confirm that extra column was not created during initial ```cube_sim_stxbp1.csv``` file processing.  If created, the extra row can be deleted manually or by adding a line in script can to remove column this automatically.



## Optional - Term Propagation
Although not necessary to run these scripts, as example CSVs are already provided, we've included a term propagation script within the [raw files directory](https://github.com/helbig-lab/STXBP1/tree/master/raw_files). This allows users to generate the base and propagated HPO terms on their own in order to view the process first hand. To run the propagation scripts:

* In order to create manual base and propogation files, change the yaml default parameters ``` pos_ic : raw_files/post_IC.csv ``` to ``` pos_ic :  ``` and then run the [R file](https://github.com/helbig-lab/STXBP1/blob/master/master_config.R) (see below) to generate a manual csv file.

```
~/Rscript master_config.R --input /path_to/input.yml
```

Note that this creates positive and negative propagation files, which are already provided in the [raw files directory](https://github.com/helbig-lab/STXBP1/tree/master/raw_files).

## Optional - Comparative Effectiveness Analysis
Although not necessary to run these scripts, we've included a set of comparative effectiveness analyses scripts, to determine the effects of certain ASMs over time. All associated raw files are included in the [raw files directory](https://github.com/helbig-lab/STXBP1/tree/master/raw_files). These allow users to generate binned information on seizure frequency and medicaiton information in order to complete this analysis:

* In order to run the comparative effectiveness analysis, ensure that the yaml file points to the correct directory for the comparative effectivness scripts. The default directly is ``` comp_dir: comp_effectiveness_analysis/ ``` 
* After confirming the directory, run the [R file](https://github.com/helbig-lab/STXBP1/blob/master/master_config.R) (see below) to generate a the comparative effectiveness results file(s).

```
~/Rscript master_config.R --input /path_to/input.yml
```

