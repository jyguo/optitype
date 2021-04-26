#!/usr/bin/env nextflow
/*

    nextflow script for running optitype on BAM files
    Copyright (c) 2014-2015 National Marrow Donor Program (NMDP)

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library;  if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

    > http://www.gnu.org/licenses/lgpl.html

  ./nextflow run nmdp-bioinformatics/flow-Optitype \
    --outfile hli-optitype.csv \
    --bamdir s3://bucket/s3/data \
    --datatype rna
*/

params.help = ''
params.datatype = 'rna'
params.optiref = "s3://gcs3-1/common/rADIO/hla_reference_rna.fasta"
outfile = file("${params.outfile}")
datatype = "${params.datatype}"
params.ncpu = 15

Channel
    .fromPath(file(params.inputManifest))
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple(row.sampleID, row.bam, row.bai) }
    .set { bam_channel }


/*  Help section (option --help in input)  */
if (params.help) {
    log.info ''
    log.info '---------------------------------------------------------------'
    log.info 'NEXTFLOW OPTITYPE'
    log.info '---------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf -with-docker nmdpbioinformatics/flow-optitype --bamdir bamfiles/ [--datatype rna] [--outfile datafile.txt]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bamdir      FOLDER             Folder containing BAM FILES'
    log.info 'Options:'
    log.info '    --datatype    STRING             Type of sequence data (default : dna)'
    log.info '    --outfile     STRING             Name of output file (default : typing_results.txt)'
    log.info ''
    exit 1
}

/* Software information */
log.info ''
log.info '---------------------------------------------------------------'
log.info 'NEXTFLOW OPTITYPE'
log.info '---------------------------------------------------------------'
log.info "Input BAM folder   (--bamdir)    : ${params.bamdir}"
log.info "Sequence data type (--datatype)  : ${params.datatype}"
log.info "Output file name   (--outfile)   : ${params.outfile}"
log.info "Project                          : $workflow.projectDir"
log.info "Git info                         : $workflow.repository - $workflow.revision [$workflow.commitId]"
log.info "\n"
log.info " parameters "
log.info " ======================"
log.info " input manifest          : ${params.inputManifest}"
log.info " output directory         : ${params.outputDir}"
log.info " ======================"
log.info ""


// Extract pair reads to fq files
process bam2fastq {
  errorStrategy 'ignore'
  tag{ subid }
  
  input:
    set val(subid), path(bamfile), path(baifile) from bam_channel
  output:
    set subid, file("${subid}.end1.fq") into fastq1
    set subid, file("${subid}.end2.fq") into fastq2

  publishDir "${params.outputDir}/${subid}"

  """
  samtools sort -@ ${params.ncpu} -n ${bamfile} | samtools bam2fq -@ ${params.ncpu} -1 ${subid}.end1.fq -2 ${subid}.end2.fq -N -O -
  """
}

//Filter the fq files
process razarEnd1 {
  errorStrategy 'ignore'
  tag{ subid }
  
  input:
    set subid, file(fq) from fastq1
    path optiref from params.optiref
  output:
    set subid, file("${subid}.raz-end1.fastq") into razarFilteredEnd1

  publishDir "${params.outputDir}/${subid}"

  """
  razers3 -tc ${params.ncpu} --percent-identity 90 --max-hits 1 --distance-range 0 --output ${subid}.raz-end1.sam ${optiref} ${subid}.end1.fq
  cat ${subid}.raz-end1.sam | grep -v ^@ | awk '{print "@"\$1"\\n"\$10"\\n+\\n"\$11}' > ${subid}.raz-end1.fastq
  """
}

//Filter the fq files
process razarEnd2 {
  errorStrategy 'ignore'
  tag{ subid }
  
  input:
    set subid, file(fq) from fastq2
    path optiref from params.optiref
  output:
    set subid, file("${subid}.raz-end2.fastq") into razarFilteredEnd2

  publishDir "${params.outputDir}/${subid}"

  """
  razers3 -tc ${params.ncpu} --percent-identity 90 --max-hits 1 --distance-range 0 --output ${subid}.raz-end2.sam ${optiref} ${subid}.end2.fq
  cat ${subid}.raz-end2.sam | grep -v ^@ | awk '{print "@"\$1"\\n"\$10"\\n+\\n"\$11}' > ${subid}.raz-end2.fastq
  """
}

//Collect filtered fq files
fqPairs = Channel.create()
fastqFiltered = razarFilteredEnd1.phase(razarFilteredEnd2).map{ fq1, fq2 -> [ fq1[0], fq1[1], fq2[1] ] }.tap(fqPairs)

//Run OptiType
process optitype {
  errorStrategy 'ignore'
  tag{ subid }

  input:
    set subid, file(fq1), file(fq2) from fastqFiltered
  output:
    file "typing_results.txt" into optioutput

  publishDir "${params.outputDir}/${subid}"

  """
  OptiTypePipeline.py -i ${fq1} ${fq2} --id ${subid} --${datatype} --outdir na > typing_results.txt
  """
}

// Print out results to output file
// optioutput
// .collectFile() {  typing ->
//       [ "typing_results.txt", typing ]
//   }
// .subscribe { file -> copy(file) }

def copy (file) { 
  log.info "Copying ${file.name} into: $outfile"
  file.copyTo(outfile)
}

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf("."))
}
