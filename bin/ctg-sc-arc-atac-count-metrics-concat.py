#!/usr/bin/python3.4

# Create multiqc-table for ARC-ATAC summary metrics

import pandas as pd
import numpy as np
import sys, getopt
import os

def main(argv):
    projdir = ''
    outputdir = ''

    usage='> Usage: ctg-sc-arc-atac-count-metrics-concat.py -i PROJECT-OUTDIR -o SUMMARY-OUTDIR'

    try:
        opts, args = getopt.getopt(argv,"hi:o:",["projdir=", "outdir="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    if len(sys.argv) <= 2:
        print("> Error: No project dir / output dir entered:")
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-i", "--projdir"):
            projdir = arg
        elif opt in ("-o", "--outdir"):
            outputdir = arg

    out1 = outputdir + "/ctg-cellranger-arc-atac-count-summary_metrics.csv"
    out2 = outputdir + "/cellranger-arc-atac-count_summary_metrics_mqc.csv"

    projid=os.path.basename(os.path.normpath(projdir))

    # list all metricsfiles
    samples = os.listdir(projdir + '/count/')

    # get id of metricfile
    sname = samples[0]
    fname = projdir + "/qc/cellranger/" + sname + ".summary.csv"
    final = pd.read_csv(fname)
    final.index = [sname]

    # concatenate all sample tables to one
    for sname in samples[1:]:
        fname = projdir + "/qc/cellranger/" + sname + ".summary.csv"
        data = pd.read_csv(fname)
        data.index = [sname]
        final = pd.concat([final,data],axis=0)

    # Write csv file
    cols = final.columns.tolist()
    cols = list( cols[i] for i in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25] )
    final = final[cols]
    final.replace(",","",regex=True)
    final.index.name = "Sample"
    final.to_csv(out1,sep=",")

    # Parse csv file to mqc
    # parse % input
    def p2f(x):
        return float(x)*100
    # parse integers with comma
    def s2i(x):
        return int(x.replace(",",""))

    mqdf = pd.read_csv(out1,
                       converters={'Estimated number of cells':s2i,
                                   'Feature linkages detected':s2i,
                                   'Linked genes':s2i,
                                   'Linked peaks':s2i,
                                   'ATAC Confidently mapped read pairs':p2f,
                                   'ATAC Fraction of genome in peaks':p2f,
                                   'ATAC Fraction of high-quality fragments in cells':p2f,
                                   'ATAC Fraction of high-quality fragments overlapping TSS':p2f,
                                   'ATAC Fraction of high-quality fragments overlapping peaks':p2f,
                                   'ATAC Fraction of transposition events in peaks in cells':p2f,
                                   'ATAC Mean raw read pairs per cell':p2f,
                                   'ATAC Median high-quality fragments per cell':p2f,
                                   'ATAC Non-nuclear read pairs':p2f,
                                   'ATAC Number of peaks':s2i,
                                   'ATAC Percent duplicates':p2f,
                                   'ATAC Q30 bases in barcode':p2f,
                                   'ATAC Q30 bases in read 1':p2f,
                                   'ATAC Q30 bases in read 2':p2f,
                                   'ATAC Q30 bases in sample index i1':p2f,
                                   'ATAC Sequenced read pairs':s2i,
                                   'ATAC TSS enrichment score':p2f,
                                   'ATAC Unmapped read pairs':p2f,
                                   'ATAC Valid barcodes':p2f
                               })

    orig_cols = mqdf.columns
    mqdf.columns = ['SampleID','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11','col12','col13','col14','col15','col16','col17','col18','col19','col20','col21','col22','col23','col24']

    f = open(out2,'a')
    f.write("# plot_type: 'table'" + "\n")
    f.write("# section_name: 'Cellranger ARC-ATAC Metrics'\n")
    f.write("# description: 'Cellranger 10x-ARC ATAC metrics'\n")
    f.write("# pconfig:\n")
    f.write("#     namespace: 'CTG ARC-ATAC'\n")
    f.write("# headers:\n")
    f.write("#     col1:\n")
    f.write("#         title: 'Sample'\n")
    f.write("#         description: 'CTG Project ID - Sample ID'\n")
    f.write("#     col2:\n")
    f.write("#         title: 'ARC Estimated Number of Cells'\n")
    f.write("#         description: 'Estimated number of cells'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#         placement: 1\n")
    f.write("#     col3:\n")
    f.write("#         title: 'ARC Feature linkages detected'\n")
    f.write("#         description: 'Feature linkages detected'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#         placement: 2\n")
    f.write("#     col4:\n")
    f.write("#         title: 'ARC Linked genes'\n")
    f.write("#         description: 'Linked genes'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#         placement: 3\n")
    f.write("#     col5:\n")
    f.write("#         title: 'ARC Linked peaks'\n")
    f.write("#         description: 'Linked peaks'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#         placement: 4\n")
    f.write("#     col6:\n")
    f.write("#         title: 'ATAC mapped (confidently) read pairs'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col7:\n")
    f.write("#         title: 'ATAC Fraction of genome in peaks'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col8:\n")
    f.write("#         title: 'ATAC Fraction of high-quality fragments in cells'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#     col9:\n")
    f.write("#         title: 'ATAC Fraction of high-quality fragments overlapping TSS'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col10:\n")
    f.write("#         title: 'ATAC Fraction of high-quality fragments overlapping peaks'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col11:\n")
    f.write("#         title: 'ATAC Fraction of transposition events in peaks in cells'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col12:\n")
    f.write("#         title: 'ATAC Mean raw read pairs per cell'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#     col13:\n")
    f.write("#         title: 'ATAC Median high-quality fragments per cell'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#     col14:\n")
    f.write("#         title: 'ATAC Non-nuclear read pairs'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col15:\n")
    f.write("#         title: 'ATAC Number of peaks'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#     col16:\n")
    f.write("#         title: 'ATAC Percent duplicates'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col17:\n")
    f.write("#         title: 'ATAC Q30 bases in barcode'\n")
    f.write("#         description: 'Q30 bases in barcode read'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col18:\n")
    f.write("#         title: 'ATAC Q30 bases in read 1'\n")
    f.write("#         description: 'Q30 bases in read 1'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col19:\n")
    f.write("#         title: 'ATAC Q30 bases in read 2'\n")
    f.write("#         description: 'Q30 bases in read 2'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col20:\n")
    f.write("#         title: 'ATAC Q30 bases in sample index i1'\n")
    f.write("#         description: 'ATAC Q30 bases in sample index i1'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col21:\n")
    f.write("#         title: 'ATAC Sequenced read pairs'\n")
    f.write("#         description: 'ATAC Sequenced read pairs'\n")
    f.write("#         format: '{:.0f}'\n")
    f.write("#     col22:\n")
    f.write("#         title: 'ATAC TSS enrichment score'\n")
    f.write("#         description: 'ATAC TSS enrichment score'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#     col23:\n")
    f.write("#         title: 'ATAC Unmapped read pairs'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col24:\n")
    f.write("#         title: 'ATAC Valid barcodes'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    mqdf.to_csv(f,sep="\t",index=False)
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
