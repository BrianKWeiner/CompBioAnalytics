#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.outDir = params.outDir ?: new Date().format('yyyy_MM_dd')

workflow {
    main:
        downloadSraFiles(params.sraList)
        convertToFastq()
        alignReads()
        countFeatures()

    // Download SRA files
    downloadSraFiles:
        Channel
            .fromFilePairs("${params.sraDir}/${params.sraList}")
            .set { sra_files }

        process download_sra {
            tag { sra_file }
            input:
                path sra_file
            output:
                path "${params.sraDir}/${params.outDir}/sra/${sra_file}"

            """
            prefetch -O ${params.sraDir}/${params.outDir} ${sra_file}
            """
        }

    // Convert SRA to FASTQ
    convertToFastq:
        sra_files
            .collectFile { file -> file.path.endsWith('.sra') }
            .set { fastq_files }

        process convert_to_fastq {
            tag { sra_file }
            input:
                path sra_file
            output:
                path "${params.sraDir}/${params.outDir}/fastq/*.fastq"

            """
            fastq-dump --outdir ${params.sraDir}/${params.outDir}/fastq --split-files ${sra_file}
            """
        }

    // Align Reads with STAR
    alignReads:
        fastq_files
            .set { bam_files }

        process run_star {
            tag { sra_file }
            input:
                path sra_file
            output:
                path "${params.sraDir}/${params.outDir}/alignments/${sra_file}Aligned.sortedByCoord.out.bam"
                path "${params.sraDir}/${params.outDir}/alignments/${sra_file}Log.out"

            """
            STAR --genomeDir ${params.genomeIndexDir} \
                 --runThreadN ${params.threads} \
                 --readFilesIn ${sra_file}_1.fastq ${sra_file}_2.fastq \
                 --outFileNamePrefix ${params.sraDir}/${params.outDir}/alignments/${sra_file} \
                 --outSAMtype BAM SortedByCoordinate \
                 --outSAMunmapped Within \
                 --outSAMattributes Standard
            """
        }

    // Count Features with featureCounts
    countFeatures:
        bam_files
            .set { count_files }

        process run_featureCounts {
            tag { sra_file }
            input:
                path bam_file
            output:
                path "${params.sraDir}/${params.outDir}/counts/${sra_file}.counts.txt"

            """
            featureCounts -a ${params.refGenomeGtf} \
                          -o ${params.sraDir}/${params.outDir}/counts/${sra_file}.counts.txt \
                          -p -T ${params.threads} \
                          ${bam_file}
            """
        }
}

