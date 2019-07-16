#!/usr/bin/env nextflow

//+++++++++++++++++ SETUP++++++++++++++++++++++++
//Name of reference organism to get RNASeq datasets for
params.reference = 'Botrytis cinerea'
params.refID = 'B.cinerea'
params.workdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/2018-05-14-Botrytis_github"
//put proteins of the reference genome here
params.refpep = "${params.workdir}/reference/*.faa.gz"
//config file for augustus
params.config = "${params.workdir}/configs/extrinsic.cfg"
//put all the scripts from the 'scripts' directory of the repository there
params.scripts = "${params.workdir}/scripts"
params.outdir = "${params.workdir}/output"
//soft masked input files for genomes we want annotated
params.input = "${params.workdir}/output/*/*soft.fasta"
//+++++++++++++++++++++++++++++++++++++++++++++++

log.info "====================================================================="
log.info "Annotation using RNASeq hints and repeat masked genomes with Augustus"
log.info "Reference : ${params.reference}"
log.info "Output    : ${params.outdir}"
log.info "====================================================================="

extrinsic_config = file(params.config)

Channel.fromPath(params.refpep)
.map{[params.refID, it]}
.tap{ compressedProteins }
.tap{ compressedProteinsForOrtho }

Channel
.fromPath(params.input)
.map{[it.simpleName, it]}
.tap{ maskedAssemblies }
.tap{ maskedAssembliesForAugustus }
.tap{ maskedAssembliesForExonerate }

process getRNASeqIDs {
    tag {params.reference}
    publishDir "${params.outdir}/", mode: 'copy'

    output:
    file 'RNAseq_ids.txt' into RNASeqIDs
    """
esearch -db taxonomy -query '${params.reference}' \
| elink -target sra \
| efetch -format docusm \
| xtract -pattern EXPERIMENT_PACKAGE \
  -if LIBRARY_STRATEGY \
  -equals 'RNA-Seq' \
  -group RUN \
  -element @accession > RNAseq_ids.txt
    """
}

RNASeqIDs
.splitText()
.map{it -> it.trim()}
.set{RNASeqIDsTrimmed}

process dumpfastq {
  publishDir "${params.outdir}/${params.refID}_fastq-dumps", mode: 'copy'
  tag { id }

  input:
    val id from RNASeqIDsTrimmed

  output:
    set id, "*.fastq" into fastqDumpForAlignment

  """
  fastq-dump $id
  """
}

process DecompressProtein1 {
  tag  { id }

  input:
  set id, "proteins.fasta.gz" from compressedProteins

  output:
  set id, "proteins.fasta" into proteinsBasic

  "zcat proteins.fasta.gz > proteins.fasta"
}

process DecompressProtein2 {
  tag  { id }

  input:
  set id, "proteins.fasta.gz" from compressedProteinsForOrtho

  output:
  set id, "proteins.fasta" into proteinsBasicForOrtho

  "zcat proteins.fasta.gz > proteins.fasta"
}

process make_exonerate_index {
  tag { id }

  input:
  set id, 'genome.fasta' from maskedAssembliesForExonerate

  output:
  set id, 'genome.fasta', 'index.esi', 'index.esd' into exn_index

  """
fasta2esd --softmask no genome.fasta index.esd
esd2esi index.esd index.esi --translate yes
  """
}

exn_index
.combine(proteinsBasic.splitFasta( by: 500, file: true))
.tap { exonerateInputs }

process exonerate {
  errorStrategy 'retry'
  maxRetries 15
  memory { 2.GB * task.attempt }
  cache 'deep'
  cpus 2
  tag { id }

  input:
  set id, 'genome.fasta', 'index.esi', 'index.esd', refID, 'altprop.fasta' from exonerateInputs

  output:
  set id, 'exn_out' into exonerate_out

  """
exonerate \
 -E false \
 --model p2g \
 --showvulgar no \
 --showalignment no \
 --showquerygff no \
 --showtargetgff yes \
 --percent 80 \
 --geneseed 250 \
 --ryo \"AveragePercentIdentity: %pi\\n\" \
 altprop.fasta \
 genome.fasta > exn_out
  """
}

// Turn the Exonerate output into hints
exonerate_out
.collectFile() { id, path -> ["${id}.hints", path.text] }
.map { path -> [path.getBaseName(), path] }
.set { exonerate_hint_input }

process exonerateOutputToHints {
  tag { id }
  publishDir "${params.outdir}/$id/hints/exonerateHints", mode: 'copy'
  input:
  set id, "exonerate.gff" from exonerate_hint_input

  output:
  set id, "out.hints" into exonerate_hints

  "${params.scripts}/exonerateGffToHints.awk < exonerate.gff > out.hints"
}

process indexAssemblyHisat2 {
  tag { id }

  input:
    set id, "genome.fasta" from maskedAssemblies

  output:
    set id, "*.ht2" into indexedAssembly

  """
  /opt/hisat2/current/hisat2-build genome.fasta $id
  """
}

process alignToAssemblyhisat2 {
  publishDir "${params.outdir}/$idAssembly/bams", mode: 'copy'
  tag { "${idAssembly} ${idFastq}" }

  input:
    set idAssembly, "${idAssembly}.*.ht2", idFastq, "${idFastq}.fastq" from indexedAssembly.combine(fastqDumpForAlignment)

  output:
    set  idAssembly, idFastq, "${idAssembly}.${idFastq}.bam" into alignedReads

  """
  /opt/hisat2/current/hisat2 -x ${idAssembly} \
  -U ${idFastq}.fastq \
  --threads 10 \
  --max-intronlen 2000 \
  | samtools view -b \
  | samtools sort \
  -o ${idAssembly}.${idFastq}.bam

  """
}

alignedReads
.tap {bamToHintsInput}
.tap {indexBamInput}

process indexBam {
  publishDir "${params.outdir}/$idAssembly/bams/", mode: 'copy'
  tag { "${idAssembly} ${idFastq}" }

  input:
  set  idAssembly, idFastq, "assembly.bam" from indexBamInput

  output:
  set  idAssembly, idFastq, "${idAssembly}.${idFastq}.bam.bai" into indexedBam

  """
  samtools index assembly.bam ${idAssembly}.${idFastq}.bam.bai
  """
}

process bamToHints {
  publishDir "${params.outdir}/$idAssembly/hints/bamHints", mode: 'copy'
  tag { "${idAssembly} ${idFastq}" }

  input:
  set idAssembly, idFastq, "assembly.bam" from bamToHintsInput

  output:
  set idAssembly, "${idAssembly}.${idFastq}.gff" into bamToHints

  """
  /opt/augustus/current/bin/bam2hints --intronsonly --maxgaplen=10 --minintronlen=15 --maxintronlen=500 --in=assembly.bam --out=${idAssembly}.${idFastq}.gff
  """
}

bamToHints
.mix(exonerate_hints)
.groupTuple(by:0)
.set{exonerateAndBamHints}

process mergeHints {
  publishDir "${params.outdir}/$idAssembly/hints", mode: 'copy'
  tag {idAssembly}

  input:
  set idAssembly, "in.*.gff" from exonerateAndBamHints

  output:
  set idAssembly, "mergedHints.gff" into mergedHints

  """
  sort -k 1,1 -k4,4n in.*.gff | awk 'BEGIN {OFS = "\t"} {\$9 = "src=E"; print \$0}' | uniq > mergedHints.gff
  """
}

mergedHints
.combine(maskedAssembliesForAugustus, by:0)
.set{mergedHintsForAugustus}

process augustusGTF {
  publishDir "${params.outdir}/$id/augustus", mode: 'copy'
  container "robsyme/augustus"
  cache 'deep'
  tag { id }

  input:
  file extrinsic_config
  set id, "hints.gff", "input.fasta" from mergedHintsForAugustus


  output:
  set id, "${id}.augustus.gff", "input.fasta" into augustusGTFoutput

  """
/opt/augustus/current/bin/augustus \
 --progress=true \
 --gff3=off\
 --softmasking=1 \
 --uniqueGeneId=true \
 --noInFrameStop=true \
 --hintsfile=hints.gff \
 --/augustus/verbosity=4 \
 --species=botrytis_cinerea \
 --extrinsicCfgFile=${extrinsic_config} \
 input.fasta > ${id}.augustus.gff
  """
}

augustusGTFoutput
.tap{augustusToFasta}
.tap{augustusToBed}


process extractFasta {
  publishDir "${params.outdir}/$id/augustus", mode: 'copy'
  tag { id }

  input:
  set id, "${id}.gff", "input.fasta" from augustusToFasta

  output:
  set id, "${id}.proteins.fasta" into augustusFastas

  """
  /opt/augustus/current/scripts/getAnnoFasta.pl --seqfile input.fasta ${id}.gff
  mv ${id}.aa ${id}.proteins.fasta
  mv ${id}.codingseq ${id}.cds.fasta
  """

}

process augustusToBed {
  publishDir "${params.outdir}/$id/augustus", mode: 'copy'
  tag { id }

  input:
  set id, "${id}.gtf", "input.fasta" from augustusToBed

  output:
  file "${id}.*.bed"

  """
  ${params.scripts}/AugustusGtfToBed.rb ${id}.gtf > ${id}.augustus.genes.bed
  ${params.scripts}/AugustusGtfToRepeatBed.rb ${id}.gtf > ${id}.augustus.repeats.bed
  ${params.scripts}/AugustusGtfToCDSpartBed.rb ${id}.gtf > ${id}.augustus.CDSpart.bed
  """
}

augustusFastas
.mix(proteinsBasicForOrtho)
.tap{ unlabledProteins }

process appendSpeciesID {
  tag { id }
  publishDir "${params.outdir}/orthoFinder/input", mode: 'copy'

  input:
  set id, "in.fasta" from unlabledProteins

  output:
  set id, "${id}.proteins.fasta" into labledProteins

  """
  awk '/^>/ {sub(/>/, ">${id}.", \$1); print \$1} /^[^>]/ {print \$0}' in.fasta > ${id}.proteins.fasta
  """
}

labledProteins
.tap{ labledProteinsForOrtho}
.tap{ labledProteinsForEffectorP }
.tap{ labledProteinsForDeepsig }
.tap{ labledProteinsForInterproscan }

process effectorP {
  tag { id }
  publishDir "${params.outdir}/$id/effectorP", mode: 'copy'

  input:
  set id, "input.fasta" from labledProteinsForEffectorP

  output:
  file "${id}.effectorP.tsv"

  """
  /opt/effectorP/current/Scripts/EffectorP.py \
  -s \
  -o ${id}.effectorP.tsv \
  -i input.fasta
  """
}

// deepsig run out of the docker container seems to have a problem using the symbolic link 'input.fasta' and crashes
// therefore I first copy it and then use the copy as input.
process deepsig {
  tag { id }
  publishDir "${params.outdir}/$id/deepsig", mode: 'copy'
  container "bolognabiocomp/deepsig"

  input:
  set id, "input.fasta" from labledProteinsForDeepsig

  output:
  file "${id}.deepsig.out"

  """
  cp input.fasta in.fasta
  docker run -v \$(pwd):/data/ bolognabiocomp/deepsig -f in.fasta -k euk -o ${id}.deepsig.out
  """

}

process interproscan {
  tag { id }
  publishDir "${params.outdir}/$id/interproscan", mode: 'copy'
  cpus 12

  input:
  set id, "proteins.fasta" from labledProteinsForInterproscan

  output:
  file "${id}.interproscan.tsv"

  """
  /opt/interproscan/current/interproscan.sh \
  --applications SignalP_EUK,Pfam,TMHMM,PANTHER,PRINTS,ProDom,ProSitePatterns,ProSiteProfiles,MobiDBLite\
  --cpu ${task.cpus} \
  --seqtype p \
  --disable-precalc \
  --goterms \
  --pathways \
  --iprlookup\
  --input proteins.fasta \
  --output-file-base ${id}.interproscan \
  --format tsv
  """

}

return

process orthoFinder {
  publishDir "${params.outdir}/orthoFinder", mode: 'copy'
  cpus 10

  input:
  file "input/*" from labledProteinsForOrtho

  output:
  file "Results_*" into orthoFinderResults

  """
  /opt/orthofinder/current/orthofinder   \
  -f input \
  -t ${task.cpus} \
  -S diamond \
  -M msa
  """
}

workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}
