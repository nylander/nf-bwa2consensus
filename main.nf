/*
- File: nf-bwa2consensus
- Last modified: 2026-04-24 12:56:16
- Sign: JN
*/

nextflow.enable.dsl=2

params.samplesheet = null
params.ref = null
params.fastq1 = null
params.fastq2 = null
params.prefix = null

workflow {

  def ch_samples

  if( params.samplesheet ) {

    ch_samples = Channel
      .fromPath(params.samplesheet)
      .splitCsv(header:true)
      .map { row ->
        def sample = (row.sample ?: "").toString()
        if( !sample ) error "samplesheet row missing 'sample'"

        def r1 = file(row.fastq1)
        def r2 = file(row.fastq2)

        def rowRef = row.containsKey('ref') ? (row.ref ?: "").toString().trim() : ""
        def refPath = rowRef ? rowRef : (params.ref ?: "")

        if( !refPath )
          error "No reference for sample '${sample}'. Provide --ref or a 'ref' column value."

        tuple(sample, r1, r2, file(refPath))
      }
  }
  else if( params.fastq1 && params.fastq2 && params.prefix ) {
    if( !params.ref ) error "Missing --ref (single-sample mode)"
    ch_samples = Channel.of(tuple(params.prefix as String, file(params.fastq1), file(params.fastq2), file(params.ref)))
  }
  else {
    error "Provide either --samplesheet OR (--fastq1 --fastq2 --prefix --ref)"
  }

  def ch_refs_keyed = ch_samples
    .map { sample, r1, r2, ref -> tuple(ref.toRealPath().toString(), ref) }
    .unique { refkey, ref -> refkey }

  def ch_ref_indexed_keyed = BWA_INDEX(ch_refs_keyed)

  def ch_samples_keyed = ch_samples
    .map { sample, r1, r2, ref -> tuple(ref.toRealPath().toString(), sample, r1, r2, ref) }

  def ch_joined = ch_samples_keyed
    .combine(ch_ref_indexed_keyed)
    .filter { sampleKey, sample, r1, r2, ref, refKey, refdir -> sampleKey == refKey }
    .map    { key, sample, r1, r2, ref, refKey, refdir -> tuple(sample, r1, r2, ref, refdir) }

  def ch_for_bwa

  if( params.fastp ) {
    def ch_reads_for_fastp = ch_joined.map { sample, r1, r2, ref, refdir -> tuple(sample, r1, r2, ref, refdir) }

    def ch_trimmed = FASTP(
      ch_reads_for_fastp.map { sample, r1, r2, ref, refdir -> tuple(sample, r1, r2) }
    ).reads

    ch_for_bwa = ch_trimmed
      .map { s, r1t, r2t -> tuple(s, tuple(r1t, r2t)) }
      .combine(ch_reads_for_fastp.map { s, r1, r2, ref, refdir -> tuple(s, tuple(ref, refdir)) })
      .filter { s1, reads, s2, refinfo -> s1 == s2 }
      .map { s, reads, s2, refinfo -> tuple(s, reads[0], reads[1], refinfo[0], refinfo[1]) }
  }
  else {
    ch_for_bwa = ch_joined
  }

  ch_bam = BWA_MEM(ch_for_bwa)

  SAMTOOLS_DEPTH(ch_bam)

  ch_bam_indexed = SAMTOOLS_INDEX_BAM(ch_bam)

  def ch_samtools_consensus = SAMTOOLS_CONSENSUS(ch_bam_indexed)
  def ch_samtools_consensus_a = SAMTOOLS_CONSENSUS_A(ch_bam_indexed)
  def ch_samtools_consensus_iupac = SAMTOOLS_CONSENSUS_IUPAC(ch_bam_indexed)

  ch_vcfgz = BCFTOOLS_MPILEUP_CALL(ch_bam_indexed)
  ch_vcfgz_indexed = BCFTOOLS_INDEX(ch_vcfgz)
  def ch_bcftools_consensus = BCFTOOLS_CONSENSUS(ch_vcfgz_indexed)

  def ch_concat_input = ch_samtools_consensus
    .map { sample_id, fasta, ref -> tuple(sample_id, fasta, ref) }
    .combine( ch_samtools_consensus_a.map { sample_id, fasta, ref -> tuple(sample_id, fasta) } )
    .filter { sid1, f1, ref, sid2, f2 -> sid1 == sid2 }
    .map { sid, f1, ref, sid2, f2 -> tuple(sid, f1, f2, ref) }
    .combine( ch_samtools_consensus_iupac.map { sample_id, fasta, ref -> tuple(sample_id, fasta) } )
    .filter { sid, f1, f2, ref, sid3, f3 -> sid == sid3 }
    .map { sid, f1, f2, ref, sid3, f3 -> tuple(sid, f1, f2, f3, ref) }
    .combine( ch_bcftools_consensus.map { sample_id, fasta -> tuple(sample_id, fasta) } )
    .filter { sid, f1, f2, f3, ref, sid4, f4 -> sid == sid4 }
    .map { sid, f1, f2, f3, ref, sid4, f4 -> tuple(sid, f1, f2, f3, f4, ref) }

  def ch_concat_out = CONCAT_CONSENSUS_FASTA(ch_concat_input)

  if( params.mafft ) {
    def ch_for_mafft = GATHER_SEQS_FOR_MAFFT(ch_concat_out)
    MAFFT_ALIGN(ch_for_mafft)
  }
}

process BWA_INDEX {
  tag { ref.baseName }

  publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, enabled: params.keepbam

  conda "bioconda::bwa=0.7.19 bioconda::samtools=1.23.1"
  container "oras://community.wave.seqera.io/library/bwa_samtools:b66f9dd49364105e"

  input:
    tuple val(refkey), path(ref)

  output:
    tuple val(refkey), path("ref_index")

  script:
  """
  set -euo pipefail
  mkdir -p ref_index
  cp -v ${ref} ref_index/${ref}
  bwa index ref_index/${ref}
  samtools faidx ref_index/${ref}
  ls -l ref_index/${ref}.*
  """
}

process BWA_MEM {
  tag { sample_id }

  publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, enabled: params.keepbam

  conda "bioconda::bwa=0.7.19 bioconda::samtools=1.23.1"
  container "oras://community.wave.seqera.io/library/bwa_samtools:b66f9dd49364105e"

  input:
  tuple val(sample_id), path(r1), path(r2), path(ref), path(ref_idx)

  output:
  tuple val(sample_id), path("${sample_id}.bam"), path(ref)

  script:
  """
  set -euo pipefail
  cp -a ${ref_idx}/* .
  bwa mem -t ${params.threads} ${ref} ${r1} ${r2} \
    | samtools sort --threads ${params.threads} -o ${sample_id}.bam
  """
}

process SAMTOOLS_INDEX_BAM {
  tag { sample_id }

  publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, enabled: params.keepbam

  conda "bioconda::samtools=1.23.1"
  container "oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a"

  input:
  tuple val(sample_id), path(bam), path(ref)

  output:
  tuple val(sample_id), path(bam), path("${bam}.bai"), path(ref)

  script:
  """
  set -euo pipefail
  samtools index -@ ${params.threads} ${bam}
  """
}

process SAMTOOLS_DEPTH {
  tag { sample_id }
  publishDir "${params.outdir}/depth", mode: 'copy', overwrite: true

  conda "bioconda::samtools=1.23.1"
  container "oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a"

  input:
  tuple val(sample_id), path(bam), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}.depth.tsv"), path(ref)

  script:
  """
  set -euo pipefail
  samtools depth -a --threads ${params.threads} ${bam} > ${sample_id}.depth.tsv
  """
}

process SAMTOOLS_CONSENSUS {
  tag { sample_id }
  publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

  conda "bioconda::samtools=1.23.1"
  container "oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a"

  input:
  tuple val(sample_id), path(bam), path(bai), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}.samtools.fasta"), path(ref)

  script:
  """
  set -euo pipefail
  samtools consensus --threads ${params.threads} ${bam} -o ${sample_id}.samtools.fasta
  tmp=${sample_id}.samtools.fasta.tmp
  sed "/^>/ s/^>\\(.*\\)/>${sample_id} samtools-consensus [ref:\\1]/" ${sample_id}.samtools.fasta > "\$tmp"
  mv "\$tmp" ${sample_id}.samtools.fasta
  """
}

process SAMTOOLS_CONSENSUS_A {
  tag { sample_id }
  publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

  conda "bioconda::samtools=1.23.1"
  container "oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a"

  input:
  tuple val(sample_id), path(bam), path(bai), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}.samtools-a.fasta"), path(ref)

  script:
  """
  set -euo pipefail
  samtools consensus --threads ${params.threads} -a ${bam} -o ${sample_id}.samtools-a.fasta
  tmp=${sample_id}.samtools-a.fasta.tmp
  sed "/^>/ s/^>\\(.*\\)/>${sample_id} samtools-a-consensus [ref:\\1]/" ${sample_id}.samtools-a.fasta > "\$tmp"
  mv "\$tmp" ${sample_id}.samtools-a.fasta
  """
}

process SAMTOOLS_CONSENSUS_IUPAC {
  tag { sample_id }
  publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

  conda "bioconda::samtools=1.23.1"
  container "oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a"

  input:
  tuple val(sample_id), path(bam), path(bai), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}.samtools-iupac.fasta"), path(ref)

  script:
  """
  set -euo pipefail
  samtools consensus --threads ${params.threads} --ambig ${bam} -o ${sample_id}.samtools-iupac.fasta
  tmp=${sample_id}.samtools-iupac.fasta.tmp
  sed "/^>/ s/^>\\(.*\\)/>${sample_id} samtools-iupac-consensus [ref:\\1]/" ${sample_id}.samtools-iupac.fasta > "\$tmp"
  mv "\$tmp" ${sample_id}.samtools-iupac.fasta
  """
}

process BCFTOOLS_MPILEUP_CALL {
  tag { sample_id }

  conda "bioconda::bcftools=1.23.1"
  container "oras://community.wave.seqera.io/library/bcftools:1.23.1--16b1a31e5dc795f7"

  input:
  tuple val(sample_id), path(bam), path(bai), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}.vcf.gz"), path(ref)

  script:
  """
  set -euo pipefail
  bcftools mpileup --fasta-ref ${ref} --annotate DP --max-depth ${params.maxdepth} ${bam} \
    | bcftools call --multiallelic-caller --output-type z \
    | bcftools view --max-alleles 2 --include 'INFO/INDEL=0 && FORMAT/DP>='${params.mindepth} \
    --output-type z --output ${sample_id}.vcf.gz
  """
}

process BCFTOOLS_INDEX {
  tag { sample_id }

  conda "bioconda::bcftools=1.23.1"
  container "oras://community.wave.seqera.io/library/bcftools:1.23.1--16b1a31e5dc795f7"

  input:
  tuple val(sample_id), path(vcfgz), path(ref)

  output:
  tuple val(sample_id), path(vcfgz), path("${vcfgz}.csi"), path(ref)

  script:
  """
  set -euo pipefail
  bcftools index --threads ${params.threads} ${vcfgz}
  """
}

process BCFTOOLS_CONSENSUS {
  tag { sample_id }
  publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

  conda "bioconda::bcftools=1.23.1"
  container "oras://community.wave.seqera.io/library/bcftools:1.23.1--16b1a31e5dc795f7"

  input:
  tuple val(sample_id), path(vcfgz), path(csi), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}.bcftools.fasta")

  script:
  """
  set -euo pipefail
  bcftools consensus --fasta-ref ${ref} --haplotype 1 --missing "N" --absent "N" ${vcfgz} \
    | awk '{print \$1}' > ${sample_id}.bcftools.fasta
  tmp=${sample_id}.bcftools.fasta.tmp
  sed "/^>/ s/^>\\(.*\\)/>${sample_id} bcftools-consensus [ref:\\1]/" ${sample_id}.bcftools.fasta > "\$tmp"
  mv "\$tmp" ${sample_id}.bcftools.fasta
  """
}

process CONCAT_CONSENSUS_FASTA {
  tag { sample_id }
  publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(samtools_fa), path(samtools_a_fa), path(samtools_iupac_fa), path(bcftools_fa), path(ref)

  output:
    tuple val(sample_id), path("${sample_id}.fasta"), path(ref)

  script:
  """
  set -euo pipefail
  cat ${samtools_fa} ${samtools_a_fa} ${samtools_iupac_fa} ${bcftools_fa} > ${sample_id}.fasta
  """
}

process GATHER_SEQS_FOR_MAFFT {
  tag { sample_id }

  input:
    tuple val(sample_id), path(consensus_fa), path(ref)

  output:
    tuple val(sample_id), path("${sample_id}.with_ref.fasta"), emit: fasta

  script:
  """
  set -euo pipefail
  cat ${ref} ${consensus_fa} > ${sample_id}.with_ref.fasta
  """
}

process MAFFT_ALIGN {
  tag { "Aligning ${sample_id}" }
  publishDir "${params.outdir}/mafft", mode: 'copy', overwrite: true

  conda 'conda-forge::mafft=7.526'
  container "community.wave.seqera.io/library/mafft:7.526--b2faf35dfc7d73ab"

  input:
    tuple val(sample_id), path(fasta)

  output:
    tuple val(sample_id), path("${sample_id}.aln.fasta"), emit: alignment

  script:
  """
  set -euo pipefail
  mafft ${params.mafft_args} ${fasta} > ${sample_id}.aln.fasta
  """
}

process FASTP {
  tag { sample_id }
  label 'process_medium'

  conda "bioconda::fastp=1.3.2"
  container "oras://community.wave.seqera.io/library/fastp:1.3.2--916946baf992e235"

  publishDir "${params.outdir}/fastp", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}_1.trimmed.fastq.gz"), path("${sample_id}_2.trimmed.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html

  script:
  """
  set -euo pipefail
  fastp \
    --thread ${task.cpus} \
    --in1 ${r1} \
    --in2 ${r2} \
    --out1 ${sample_id}_1.trimmed.fastq.gz \
    --out2 ${sample_id}_2.trimmed.fastq.gz \
    --json ${sample_id}.fastp.json \
    --html ${sample_id}.fastp.html \
    --detect_adapter_for_pe \
    ${params.fastp_args}
  """
}
