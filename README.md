# ctg-sc-arc-10x
## nextflow pipeline for 10x-sc-ARC (10x multiome RNA+ATAC) analysis with cellranger-arc

- Analyze 10x arc-seq (RNA + ATAC) in one pipeline. 
- Demux -> cellranger-arc -> QC
- Supports ATAC and RNA libraries sequenced on same flowcell or two different flowcells.
- Supports different indexing of RNA and ATAC library (e.g. RNA dual and ATAC single). See `Handle dual and single indexing in same sequencing run` for more info.


1. Clone and build the Singularity container for this pipeline: https://github.com/perllb/ctg-sc-arc-10x/tree/master/container
2. Edit your samplesheet to match the example samplesheet. See section `SampleSheet` below
3. Edit the nextflow.config file to fit your project and system. 
4. Run pipeline 
```
nohup nextflow run pipe-sc-arc-10x.nf > log.pipe-sc-arc-10x.txt &
```

## Input files

1. Samplesheet (CTG_SampleSheet.sc-arc-10x.csv)

### Samplesheet requirements:

Note: no header! only the rows shown below, starting with the column names.
Note: Must be in comma-separated values format (.csv)

 | Sample_ID | index | Sample_Project | Sample_Species | Sample_Lib | Sample_Pair | 
 | --- | --- | --- | --- | --- | --- | 
 | Sr1 | SI-GA-D9 | proj_2021_012 | human | rna | 1 |
 | Sr2 | SI-GA-H9 | proj_2021_012 | human | rna | 2 |
 | Sat1 | SI-GA-C9 | proj_2021_012 | human | atac | 1 |
 | Sat2 | SI-GA-C9 | proj_2021_012 | human | atac | 2 |

- The nf-pipeline takes the following Columns from samplesheet to use in channels:


- `Sample_ID` : ID of sample. Sample_ID can only contain a-z, A-Z and "_".  E.g space and hyphen ("-") are not allowed! If 'Sample_Name' is present, it will be ignored. 
- `index` : Must use index ID (10x ID) if dual index. For single index, the index sequence works too.
- `Sample_Project` : Project ID. E.g. 2021_033, 2021_192.
- `Sample_Species` : Only 'human'/'mouse'/'custom' are accepted. If species is not human or mouse, set 'custom'. This custom reference genome has to be specified in the nextflow config file. See below how to edit the config file.
- `Sample_Lib` : 'rna'/'atac'. Specify whether sample is rna or atac library. 
- `Sample_Pair` : To match the rna sample with the corresponding atac sample. e.g. in the example above, sample 'Sr1' is the rna library, that should be matched with 'Sat1' which is the atac library of the sample

### CSV format templates

#### 1. Samplesheet : `CTG_SampleSheet.sc-arc-10x.csv`
```
Sample_ID,index,Sample_Project,Sample_Species,Sample_Lib,Sample_Pair
Si1,SI-GA-D9,2021_012,human,rna,1
Si2,SI-GA-H9,2021_012,mouse,rna,2
Sample1,SI-GA-C9,2021_013,human,atac,1
Sample2,SI-GA-C9,2021_013,mouse,atac,2
``` 

## Pipeline steps:

Cellranger-arc version: cellranger-arc v2.0.0

* `parse samplesheets`: Creates samplesheets (one for RNA, and one for ATAC) for demux based on the input samplesheet. 
* `generate library csv`: Creates library.csv needed for `count`. Based on input samplesheet. One .csv per matched RNA and ATC sample.
* `Demultiplexing` (cellranger-arc mkfastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/mkfastq). Does this separately for RNA and ATAC (since they often have different index types (dual/single) or are sequenced on different flow cells)
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `Align` + `Counts` (cellranger-arc count): Aligns fastq files to reference genome, counts genes for each cell/barcode, and perform ATAC analysis - Then performs secondary analysis such as clustering and generates the cloupe files (https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/count).
* `Aggregation` (cellranger-arc aggr): Automatically creates the input csv pointing to atac_fragments, per_barcode_metrics, gex_molecule_info files for each sample to be aggregated and executes aggregation (https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/aggr). 
* `Cellranger count metrics` (bin/ctg-sc-arc-atac-seq-count-metrics-concat.py and bin/ctg-sc-arc-gex-seq-count-metrics-concat.py): Collects main count metrics from gex and atac libraries (#cells, #reads/cell, #peaks  etc.) from each sample and collect in table that is presented by multiqc.
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report
* `md5sum`: md5sum of all generated files



## Handle dual and single indexing in same sequencing run

If your RNA and ATAC libraries have different indexing it can be handled as following:

#### RNA dual - ATAC single
In nextflow.config, set 
```
// bcl2fastq arguments
bcl2fastqarg_rna = "" 
bcl2fastqarg_atac = "--use-bases-mask=Y28n*,I6n*,N10,Y90n*" 
// Index type ('dual' or 'single')
index_rna = "dual"
index_atac = "single"	

```
	
## Container
- `sc-arc-seq-10x`: For 10x sc-arc-seq. Based on cellranger-arc v2.0.0.
https://github.com/perllb/ctg-sc-arc-10x/tree/master/container

Build container:
NOTE: Environment.yml file has to be in current working directory
```
sudo -E singularity build sc-arc-seq-10x.sif sc-cite-seq-10x-builder 
```

Add path to .sif in nextflow.config

## Output:
* ctg-PROJ_ID-output
    * `qc`: Quality control output. 
        * cellranger metrics: Main metrics summarising the count / cell output 
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * multiqc output: Summarizing FastQC output and demultiplexing (https://multiqc.info/)
    * `fastq`: Contains raw fastq files from cellranger-arc mkfastq.
    * `count`: Cellranger-arc count output. Here you find gene/cell count matrices, atac_fragments, peaks, feature quantification, secondary analysis output, and more. See (https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/count) for more information on the output files.
    * `summaries`: 
        * web-summary files which provide an overview of essential metrics from the 10x run. 
        * cloupe files which can be used to explore the data interactively in the Loupe browser (https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/visualization/latest/what-is-loupe-browser)  
    * `aggregate`:
        * Output from cellranger aggregation. 
    * `ctg-md5.PROJ_ID.txt`: text file with md5sum recursively from output dir root    



## Custom genome 

If custom genome (not hg38 or mm10) is used

1. Set "Sample_Species" column to 'custom' in samplesheet:

Example:
 | Sample_ID | index | Sample_Project | Sample_Species | Sample_Lib | Sample_Pair | 
 | --- | --- | --- | --- | --- | --- | 
 | Sr1 | SI-GA-D9 | proj_2021_012 | **custom** | rna | 1 |
 | Sr2 | SI-GA-H9 | proj_2021_012 | **custom** | rna | 2 |
 | Sat1 | SI-GA-C9 | proj_2021_012 | **custom** | atac | 1 |
 | Sat2 | SI-GA-C9 | proj_2021_012 | **custom** | atac | 2 |
 
 2. In nextflow.config, set 
 `custom_genome=/PATH/TO/CUSTOMGENOME`
 
## Add custom genes (e.g. reporters) to cellranger annotation

Use the `ctg-cellranger-add2ref` script. 

https://github.com/perllb/ctg-cellranger-add2ref

