# RNA-seq Analysis pipeline, and Differential Expression and Downstream Analysis (NGS)

This project explores the reasons for wastage of potatoes (*Solanum tuberosum*) caused by post-harvest related to breaking dormancy and premature sprouting on a molecular level. 

**Data:** The data contains 6 files of RNA-seq raw sequence Paired-End reads, under two conditions; (Week 0 - at postharvest; and at Week 1  of storage and just before the sprouting - Three replicates for each and all reads corresponding to chromosome 12). As follows:


|  Week & Replicate        |  R1 Reads             |  R1 Reads             |
| ------------------------ | --------------------- | --------------------- |
| **Week 0- Replicate 1:** | ST_w0_R1_chr12.r1.fq  | ST_w0_R1_chr12.r2.fq  |
| **Week 0- Replicate 2:** | ST_w0_R2_chr12.r1.fq  | ST_w0_R2_chr12.r2.fq  |
| **Week 0- Replicate 3:** | ST_w0_R3_chr12.r1.fq  | ST_w0_R3_chr12.r2.fq  |
| **Week 1- Replicate 1:** | ST_w0_R1_chr12.r1.fq  | ST_w0_R1_chr12.r2.fq  |
| **Week 1- Replicate 2:** | ST_w0_R2_chr12.r1.fq  | ST_w0_R2_chr12.r2.fq  |
| **Week 1- Replicate 3:** | ST_w0_R3_chr12.r1.fq  | ST_w0_R3_chr12.r2.fq  |

The data is added to the repository.

**Objectives:** The objective is to develop a bioinformatics RNA-seq analysis pipeline using best practices, and perforrm differential expression and downstream analysis. 


## Requirements

To run the Linux PBS job you will need the latest versions of:
* `FastQC`
* `STAR` aligner
* `HTSeq`

To run the DIfferential expression and downstream analysis, you will need:

* Latest version of `R`. 
* Any relevant function are attached within this repository.


## The Pipeline
The pipeline split into two parts; 1) Bash PBS job, 2) Differential expression and downstream analysin using R.

### Part I

Bash script (`Complete_Script.sub`) to automate the Quality control, read alignment and conting the number of reads map to each feature. This script can be submitted as a PBS job submission as it is computaiotnally intensive.

**Note: User must replace the `/PATH/TO/` in line 19 (`tar xvf /PATH/TO/dataForAssignment.tar.gz`) with the absolute path of the directory where the data stored.**

The automated bash PBS script performs tasks below once submitted:
1. Makes a new directory within the working directory.
2. Copies the required data into newly created directory.
3. Changes into the directory with data once the data transfer is complete.
4. Creates a new directory to store the `FastQC_Reports`.
5. changes into `FastQC_Reports` directory.
5. Loads the `FastQC` module.
6. Generates the `FastQC` reports and saves them in the parent directory of the current directory.
7. Purges the `FastQC` module and changes into its parent directory.
8. Loads the `STAR` aligner module.
9. Creates 6 `STAR_Output_*` directories to store the alignment output.
10. Creates a new directory for `STAR` index.
11. Runs the `STAR` aligner using the parameters provided.
12. Purges the `STAR` module.
13. Loads the `HTSeq` modules.
14. Runs the `HTSeq` using parameters provided.
15. Purges the `HTSeq` modulw.
16. Exit the PBS job.

The output form bash pipeline is used in differential expression and downstream analysis.


### Part II

In this section, `R` is used to perform differential expression (DE) between two conditions and downstream analysis. The main R program is saved as `limma_DEAnalysis.R` and uses `limma` package. Output files (annotation, KEGG pathways) will be saved in working directory. This program performs tasks listed below:

1. Once the data and required packages are loaded, it generates DE plot with "Voom". 
2. A multivaariate analysis analysis showing how good/bad the replicates within each conditionare clustered.
    * PCA
    * HCA
3. Reports the "total number" of over- and under-expressed genes.
4. Produces an annotated list of the "top 100" DE genes and saves it as  `potato_annot.csv` file.
5. Produces an "enhanced volcano" plot with top 10 DE genes highlighted.
6. Fetches the KEGG metabolic pathways for top 10 DE genes and integrates them into the gene expression list. The output is      saveed as `kegg_pathways_results.csv`.







