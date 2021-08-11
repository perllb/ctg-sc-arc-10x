#!/usr/bin/env nextFlow

// Base params
runfolder_rna = params.runfolder_rna
runfolder_atac = params.runfolder_atac
basedir = params.basedir
metaid = params.metaid

// Output dirs
outdir = params.outdir
fqdir = params.fqdir
qcdir = params.qcdir
countdir = params.countdir
aggdir = params.aggdir
ctgqc = params.ctgqc
metadir = params.metadir

// Demux args
b2farg_rna = params.bcl2fastqarg_rna
b2farg_atac = params.bcl2fastqarg_atac
index_rna = params.index_rna
index_atac = params.index_atac
demux = params.demux

// Read and process CTG samplesheet 
sheet = file(params.sheet)

// create new samplesheet in cellranger mkfastq IEM (--samplesheet) format. This will be used only for demultiplexing
newsheet_rna = "$metadir/samplesheet.nf.sc-arc-10x.rna.csv"
newsheet_atac = "$metadir/samplesheet.nf.sc-arc-10x.atac.csv"

println "============================="
println ">>> sc-arc-10x pipeline >>>"
println ""
println "> INPUT: "
println ""
println "> run-meta-id		: $metaid "
println "> basedir		: $basedir "
println "> runfolder rna	: $runfolder_rna "
println "> runfolder atac	: $runfolder_atac "
println "> sample-sheet		: $sheet "
println ""
println " - demux settings " 
println "> bcl2fastq-arg-rna    : '${b2farg_rna}' "
println "> bcl2fastq-arg-atac   : '${b2farg_atac}' "
println "> demux                : $demux " 
println "> index-rna            : $index_rna "
println "> index-atac           : $index_atac "
println "> metadir              : $metadir "
println ""
println " - output directory structure "
println "> outdir               : $outdir "
println "> fastq                : $fqdir "
println "> qc                   : $qcdir "
println "> count                : $countdir " 
println "> aggregated           : $aggdir "
println "> ctg-qc               : $ctgqc "
println "> metadata             : $metadir "
println ""
println "============================="


// extract RNA samplesheet info
Channel
    .fromPath(sheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species, row.Sample_Lib, row.Sample_Pair ) }
    .tap{infoall}
    .into { crlib_ch; cragg_ch; fqc_ch }

// RNA Projects
Channel
    .fromPath(sheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .set { count_summarize  }

println " > Samples to process: "
println "[Sample_ID,Sample_Project,Sample_Species,Sample_Lib,pair]"
infoall.subscribe { println "Info: $it" }

println " > RNA Projects to process : "
println "[Sample_Project]"
infoProject.subscribe { println "Info Projects: $it" }

// Parse RNA samplesheet 
process parsesheet_rna {

	tag "$metaid"

	input:
	val sheet
	val index_rna

	output:
	val newsheet_rna into demux_sheet_rna

	when:
	demux == 'y'

	"""
mkdir -p $metadir
cat $sheet | grep \',rna,\\|Lane,Sample_ID\' > tmp.sheet.RNA.csv
python $basedir/bin/ctg-parse-samplesheet.10x.py -s tmp.sheet.RNA.csv -o $newsheet_rna -i $index_rna
rm tmp.sheet.RNA.csv
	"""
}

// Parse ATAC samplesheet 
process parsesheet_atac {

	tag "$metaid"

	input:
	val sheet
	val index_atac

	output:
	val newsheet_atac into demux_sheet_atac

	when:
	demux == 'y'

	"""
mkdir -p $metadir
cat $sheet | grep \',atac\\|Lane,Sample_ID\' > tmp.sheet.ATAC.csv
python $basedir/bin/ctg-parse-samplesheet.10x.py -s tmp.sheet.ATAC.csv -o $newsheet_atac -i $index_atac
rm tmp.sheet.ATAC.csv
	"""
}

// aggregation
process gen_libraries_csv {

    	tag "${sid}_${projid}"

	input:
	val sheet
	set sid, projid, ref, lib, pair from crlib_ch

 	output:
	set sid, projid, ref, lib, pair into count_lib_csv

	when:
	lib == 'rna'

	"""
mkdir -p $metadir

libcsv=$metadir/${projid}_${sid}_libraries.csv

# Print header
echo 'fastqs,sample,library_type' > \$libcsv
# Print RNA entry
echo '${fqdir}/rna/${projid},$sid,Gene Expression' >> \$libcsv
# Get paired ATAC sample
atacid=\$(grep ',atac,$pair' $sheet | cut -f2 -d ',')
echo "${fqdir}/atac/${projid},\$atacid,Chromatin Accessibility" >> \$libcsv

        """
}

// Run RNA mkFastq
process mkfastq_rna {

	tag "${metaid}-rna"

	input:
	val rnasheet from demux_sheet_rna
	
	output:
	val "count" into count_rna
	val "gorna" into fastqc_go_rna

	when:
	demux == 'y'

	"""
if [ '$index_rna' == 'dual' ]; then
   indexarg='--filter-dual-index'
else
   indexarg='--filter-single-index'
fi

cellranger-arc mkfastq \\
	   --id=${metaid}_rna \\
	   --run=$runfolder_rna \\
	   --samplesheet=$rnasheet \\
	   --jobmode=local \\
	   --localmem=150 \\
	   --output-dir ${fqdir}/rna \\
	   --delete-undetermined \\
	   \${indexarg} \\
	   $b2farg_rna

"""
}

// Run ATAC mkFastq
process mkfastq_atac {

	tag "${metaid}-atac"

	input:
	val atacsheet from demux_sheet_atac

	output:
	val "count" into count_atac
	val "goatac" into fastqc_go_atac

	when:
	demux == 'y'

	"""
if [ '$index_atac' == 'dual' ]; then
   indexarg='--filter-dual-index'
else
   indexarg='--filter-single-index'
fi

cellranger-arc mkfastq \\
	   --id=${metaid}_atac \\
	   --run=$runfolder_atac \\
	   --samplesheet=$atacsheet \\
	   --jobmode=local \\
	   --localmem=150 \\
	   --output-dir ${fqdir}/atac \\
	   --delete-undetermined \\
	   \${indexarg} \\
	   $b2farg_atac
"""
}

// count RNA + ATAC
process count {

	tag "${sid}-${projid}"
	publishDir "${countdir}/", mode: "copy", overwrite: true

	input: 
	val rnaready from count_rna
	val atacready from count_atac
	set sid, projid, ref, lib, pair from count_lib_csv

	output:
	file "${sid}/outs/" into samplename
	val "${qcdir}/cellranger/${sid}.summary.csv" into count_metrics
	val "${aggdir}/${sid}.molecule_info.h5" into count_agg

	when:
	lib == 'rna'

	"""
if [ $ref == "Human" ] || [ $ref == "human" ]
then
	genome="/projects/fs1/shared/references/hg38/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
	genome="/projects/fs1/shared/references/mm10/cellranger/refdata-cellranger-arc-mm10-2020-A-2.0.0"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ] 
then
	genome=${params.custom_genome}
else
	echo ">SPECIES NOT RECOGNIZED!"
	genome="ERR"
fi

mkdir -p ${countdir}

libcsv=$metadir/${projid}_${sid}_libraries.csv

cellranger-arc count \\
	--id=$sid \\
	--libraries=\$libcsv \\
	--reference=\$genome \\
	--localmem=150 \\
	--jobmode=local \\
	--localcores=${task.cpus} 


## Copy files for aggregation
# h5 file 
mkdir -p $aggdir
cp ${sid}/outs/gex_molecule_info.h5 ${aggdir}/${sid}.gex_molecule_info.h5

# Atac fragments
cp ${sid}/outs/atac_fragments.tsv.gz ${aggdir}/${sid}.atac_fragments.tsv.gz
cp ${sid}/outs/atac_fragments.tsv.gz.tbi ${aggdir}/${sid}.atac_fragments.tsv.gz.tbi

# Per Barcode metrics
cp ${sid}/outs/per_barcode_metrics.csv ${aggdir}/${sid}.per_barcode_metrics.csv 

## Copy metrics file for qc
# Remove if it exists
if [ -f ${qcdir}/cellranger/${sid}.summary.csv ]; then
	rm -r ${qcdir}/cellranger/${sid}.summary.csv
fi
mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger/
cp ${sid}/outs/summary.csv ${qcdir}/cellranger/${sid}.summary.csv

## Copy to delivery folder 
mkdir -p ${outdir}/summaries
mkdir -p ${outdir}/summaries/cloupe
mkdir -p ${outdir}/summaries/web-summaries
cp ${sid}/outs/web_summary.html ${outdir}/summaries/web-summaries/${sid}.web_summary.html
cp ${sid}/outs/cloupe.cloupe ${outdir}/summaries/cloupe/${sid}_cloupe.cloupe

## Copy to CTG-QC dir 
mkdir -p ${ctgqc}
mkdir -p ${ctgqc}/web-summaries
cp ${sid}/outs/web_summary.html ${ctgqc}/web-summaries/${sid}.web_summary.html

	"""

}

process summarize_count {

	tag "${projid}"

	input:
	val metrics from count_metrics.collect()

	output:
	val "y" into mqc_count 	
	val "x" into run_summarize

	"""
cd $outdir

mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger

# RNA (GEX) summaries
#python $basedir/bin/ctg-sc-arc-gex-count-metrics-concat.py -i ${outdir} -o ${qcdir}/cellranger
# ATAC summaries
#python $basedir/bin/ctg-sc-arc-atac-count-metrics-concat.py -i ${outdir} -o ${qcdir}/cellranger
	"""
}

// aggregation
process gen_aggCSV {

	tag "${sid}_${projid}"

	input:
	set sid, projid, ref, lib, pair from cragg_ch

	output:
	set projid, ref into craggregate

	when:
	lib == 'rna'

	"""
mkdir -p ${aggdir}
aggcsv=${aggdir}/${projid}_libraries.csv

if [ -f \${aggcsv} ]
then
	if grep -q $sid \$aggcsv
	then
		echo ""
	else
		echo "${sid},${aggdir}/${sid}.atac_fragments.tsv.gz,${aggdir}/${sid}.per_barcode_metrics.csv,${aggdir}/${sid}.gex_molecule_info.h5" >> \$aggcsv
	fi
else
	echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > \$aggcsv
	echo "${sid},${aggdir}/${sid}.atac_fragments.tsv.gz,${aggdir}/${sid}.per_barcode_metrics.csv,${aggdir}/${sid}.gex_molecule_info.h5" >> \$aggcsv

fi


	"""
}

process aggregate {

	publishDir "${outdir}/aggregate/", mode: 'move', overwrite: true
	tag "$projid"
  
	input:
	set projid, ref from craggregate.unique()
	val moleculeinfo from count_agg.collect()

	output:
	file "${projid}_agg/outs" into doneagg
	val projid into md5_agg_go

	"""
if [ $ref == "Human" ] || [ $ref == "human" ]
then
	genome="/projects/fs1/shared/references/hg38/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
	genome="/projects/fs1/shared/references/mm10/cellranger/refdata-cellranger-arc-mm10-2020-A-2.0.0"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ] 
then
	genome=${params.custom_genome}
else
	echo ">SPECIES NOT RECOGNIZED!"
	genome="ERR"
fi


cellranger-arc aggr \
   --id=${projid}_agg \
   --csv=${aggdir}/${projid}_libraries.csv \
   --reference=\$genome

## Copy to delivery folder 
cp ${projid}_agg/outs/web_summary.html ${outdir}/summaries/web-summaries/${projid}_agg.web_summary.html
mv ${projid}_agg/outs/cloupe.cloupe ${outdir}/summaries/cloupe/${projid}_agg_cloupe.cloupe

## Copy to CTG QC dir 
cp ${outdir}/summaries/web-summaries/${projid}_agg.web_summary.html ${ctgqc}/web-summaries/

## Remove the gex_molecule_info.h5 files that are stored in the aggregate folder (the original files are still in count-cr/../outs 
rm ${aggdir}/*h5

## Remove the barcode_metrics.csv files that are stored in the aggregate folder (the original files are still in count-cr/../outs 
rm ${aggdir}/*barcode_metrics.csv

## Remove the atac_fragments.csv files that are stored in the aggregate folder (the original files are still in count-cr/../outs 
rm ${aggdir}/*atac_fragments.tsv.gz*

	"""
}

process fastqc {

	tag "${sid}-${projid}"

	input:
	val gorna from fastqc_go_rna
	val goatac from fastqc_go_atac
	set sid, projid, ref, lib, pair from fqc_ch	
		
	output:
	val projid into mqc_cha
	val projid into md5_fastqc_go

	"""
mkdir -p ${qcdir}
mkdir -p ${qcdir}/fastqc

for file in ${fqdir}/${lib}/${projid}/${sid}*fastq.gz
	do fastqc -t ${task.cpus} \$file --outdir=${qcdir}/fastqc
done
	"""
}

process multiqc_count_run {

	tag "${metaid}"

	input:
	val x from run_summarize.collect()
	val projid from mqc_cha.collect()
		
	output:
	val "x" into summarized

	"""
# Edit atac stats.json to fetch stats in multiqc
# Get flowcell id
cd $fqdir/atac/Stats
flow=\$(grep Flowcell Stats.json | cut -f4 -d\"\\"\")
sed "s/\$flow/\${flow}_ATAC/g" Stats.json > tmp.txt
mv tmp.txt Stats.json 

cd ${outdir}
mkdir -p ${qcdir}/multiqc
multiqc -f ${fqdir} ${qcdir}/fastqc/ ${qcdir}/cellranger/ ${fqdir}/atac/Stats/Stats.json ${fqdir}/rna/Stats/Stats.json --outdir ${qcdir}/multiqc -n ${metaid}_sc-aft-rna-10x_summary_multiqc_report.html

cp -r ${qcdir} ${ctgqc}

	"""

}

// Final process, when all is done: md5 recursively from output root folder
process md5sum {

	tag "${projid}"

	input:
	val v from summarized.collect()
	val projid from md5_fastqc_go.unique()
	val x from md5_agg_go.collect()
	
	"""
# Remove Undetermined files!
rm ${fqdir}/atac/Undetermined*
rm ${fqdir}/rna/Undetermined*

cd ${outdir}
find -type f -exec md5sum '{}' \\; > ctg-md5.${projid}.txt

touch $runfolder/ctg.sc-arc-10x.done
        """ 

}
