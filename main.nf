#!/usr/bin/env nextflow

/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 'nf-run-snpAD' - A Nextflow pipeline to run snpAD
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Judith Ballesteros Villascán
 GitHub: https://github.com/jbv2/run-snpAD
 ----------------------------------------------------------------------------------------
 */

/* 
 Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */ 

params.input_tsv = null  // This will hold the TSV file path from the command line

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BAM2SNPAD_NOT_UDG } from './modules/local/bam2snpad/not_udg/main'
include { BAM2SNPAD_UDG     } from './modules/local/bam2snpad/udg/main'
include { MAPABILITY        } from './modules/local/mapability/main'
include { ACCESIBILITY      } from './modules/local/accesibility/main'
include { SNPAD             } from './modules/local/snpAD/main'
include { CALL              } from './modules/local/call/main'
include { BCFTOOLS_CONCAT   } from './modules/local/bcftools/concat/main'
include { BCFTOOLS_FILTER   } from './modules/local/bcftools/filter/main'
include { BCFTOOLS_SORT     } from './modules/local/bcftools/sort/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SORTED    } from './modules/local/bcftools/index/main'
include { BCFTOOLS_ANNOTATE     } from './modules/local/bcftools/annotate/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_ANNOTATED    } from './modules/local/bcftools/index/main'
include { BCFTOOLS_STATS   } from './modules/local/bcftools/stats/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


workflow {

    if (!params.input_tsv) {
        error "Please provide a TSV file using the '--input_tsv' parameter."
    }

    // Read the TSV file and extract columns
    ch_samples = Channel
        .fromPath(params.input_tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def meta = [ id: row.sample_id ]
            def bam  = row.inputbam
            def bai  = row.bai
            def udg  = row.udg
            return [meta, bam, bai, udg]
        }
    
    // Filter samples where UDG is 'true'
    ch_udg_samples = ch_samples.filter { meta, bam, bai, udg -> udg == "true" }
    ch_non_udg_samples = ch_samples.filter { meta, bam, bai, udg -> udg == "false" }

    // Separate UDG and non-UDG samples
    ch_ref_fasta = Channel.from(params.ref_fasta)
    ch_i = Channel.of(1..22, 'X', 'Y')

    ch_input_bam2snpad_udg = ch_udg_samples
    .combine(ch_ref_fasta)
    .combine(ch_i)
    .multiMap{
            meta, bam, bai, udg, fasta, chrom ->
            bam:      [ meta, bam ]
            bai:      bai       
            fasta:    fasta     
            chrom:    chrom 
        }

    ch_input_bam2snpad_non_udg = ch_non_udg_samples
    .combine(ch_ref_fasta)
    .combine(ch_i)
    .multiMap{
            meta, bam, bai, udg, fasta, chrom ->
            bam:      [ meta, bam ]
            bai:      bai       
            fasta:    fasta     
            chrom:    chrom 
        }
    
    ch_input_not_udg = BAM2SNPAD_NOT_UDG(ch_input_bam2snpad_non_udg.bam, ch_input_bam2snpad_non_udg.bai, ch_input_bam2snpad_non_udg.fasta, ch_input_bam2snpad_non_udg.chrom.flatten())
    ch_input_udg = BAM2SNPAD_UDG(ch_input_bam2snpad_udg.bam, ch_input_bam2snpad_udg.bai, ch_input_bam2snpad_udg.fasta, ch_input_bam2snpad_udg.chrom.flatten())

    // Combine UDG and non-UDG processed samples
    ch_input = Channel.empty()
    if (ch_input_not_udg) {
    ch_input = ch_input.mix(ch_input_not_udg)
    }
    if (ch_input_udg) {
    ch_input = ch_input.mix(ch_input_udg)
    }

    ch_map_bed = Channel.from(params.map_bed) 

    ch_input_mappability = ch_input
    .combine(ch_map_bed)
    .multiMap{
            meta, input, chrom, map ->
            input:    [ meta, input ]
            chrom:    [ chrom       ]
            map_bed:  [ map         ]
        }

    ch_mapped_inputs = MAPABILITY(ch_input_mappability.input, ch_input_mappability.chrom.flatten(), ch_input_mappability.map_bed)

    // Add accesibility masked regions 1KGP

    if ( params.accesibility ) {
    //Run Accesibility
    ch_accesibility_bed = Channel.from(params.accesibility_bed)

    ch_input_accesibility = ch_mapped_inputs
    .combine(ch_accesibility_bed)
    .multiMap{
            meta, input, chrom, acces_bed ->
            input:     [ meta, input  ]
            chrom:     [ chrom        ]
            acces_bed: [ acces_bed    ]
        }
    
    ch_accesibility_inputs = ACCESIBILITY(ch_input_accesibility.input, ch_input_accesibility.chrom.flatten(), ch_input_accesibility.acces_bed)

    // Continue with the same further steps
    ch_snpad_input = ch_accesibility_inputs
    .filter { entry ->
    entry[2] == 21 //Here use chrom 21
    }
    .map { meta, input, chrom ->
    [ meta, input ]
    }
    
    ch_snpad_params =  SNPAD(ch_snpad_input)

    ch_ref_fasta_fai = Channel.from(params.ref_fasta_fai)

    ch_call_input = SNPAD.out.priors
    .combine(SNPAD.out.errors, by: 0)
    .combine(ch_accesibility_inputs, by: 0) 
    .combine(ch_ref_fasta_fai)
    .multiMap {
        meta, priors, errors, input, chrom, fai ->
        priors:   [ priors      ]
        errors:   [ errors      ]
        input:    [ meta, input ]
        chrom:    [ chrom       ]
        fai:      [ fai         ]
            }

    ch_split_vcfs = CALL(ch_call_input.priors, ch_call_input.errors, ch_call_input.input, ch_call_input.chrom.flatten(), ch_call_input.fai.flatten())
    .groupTuple(size : 24)
    // .combine(ch_ref_fasta)
    .multiMap{ meta, vcfs ->
        vcf:   [meta, vcfs]
        }

    ch_concatenated = BCFTOOLS_CONCAT(ch_split_vcfs)
    .combine(ch_ref_fasta)
    .multiMap{ meta, vcfs, fasta ->
        vcf:   [meta, vcfs]
        fasta: [fasta]
        }

    ch_filtered = BCFTOOLS_FILTER(ch_concatenated)

    ch_sorted = BCFTOOLS_SORT(ch_filtered) 

    ch_index_sorted = BCFTOOLS_INDEX_SORTED(ch_sorted)
    
    ch_input_annotate = ch_sorted
    .combine(ch_index_sorted, by: 0)
    .multiMap{ meta, vcf, tbi ->
    vcfs: [meta, vcf]
    tbi:  [tbi]
    }

    ch_annotate = BCFTOOLS_ANNOTATE(ch_input_annotate)

    ch_annotate_index = BCFTOOLS_INDEX_ANNOTATED(ch_annotate)

    ch_stats_input = ch_annotate
    .combine(ch_annotate_index, by: 0)

    BCFTOOLS_STATS(ch_stats_input)
    
    } else {
    // Continue without accesibility 
    ch_snpad_input = ch_mapped_inputs
    .filter { entry ->
    entry[2] == 21 //Here use chrom 21
    }
    .map { meta, input, chrom ->
    [ meta, input ]
    }
    
    ch_snpad_params =  SNPAD(ch_snpad_input)

    ch_ref_fasta_fai = Channel.from(params.ref_fasta_fai)

    ch_call_input = SNPAD.out.priors
    .combine(SNPAD.out.errors, by: 0)
    .combine(ch_mapped_inputs, by: 0) 
    .combine(ch_ref_fasta_fai)
    .multiMap {
        meta, priors, errors, input, chrom, fai ->
        priors:   [ priors      ]
        errors:   [ errors      ]
        input:    [ meta, input ]
        chrom:    [ chrom       ]
        fai:      [ fai         ]
            }

    ch_split_vcfs = CALL(ch_call_input.priors, ch_call_input.errors, ch_call_input.input, ch_call_input.chrom.flatten(), ch_call_input.fai.flatten())
    .groupTuple(size : 24)
    // .combine(ch_ref_fasta)
    .multiMap{ meta, vcfs ->
        vcf:   [meta, vcfs]
        }

    ch_concatenated = BCFTOOLS_CONCAT(ch_split_vcfs)
    .combine(ch_ref_fasta)
    .multiMap{ meta, vcfs, fasta ->
        vcf:   [meta, vcfs]
        fasta: [fasta]
        }

    ch_filtered = BCFTOOLS_FILTER(ch_concatenated)

    ch_sorted = BCFTOOLS_SORT(ch_filtered) 

    ch_index_sorted = BCFTOOLS_INDEX_SORTED(ch_sorted)
    
    ch_input_annotate = ch_sorted
    .combine(ch_index_sorted, by: 0)
    .multiMap{ meta, vcf, tbi ->
    vcfs: [meta, vcf]
    tbi:  [tbi]
    }

    ch_annotate = BCFTOOLS_ANNOTATE(ch_input_annotate)

    ch_annotate_index = BCFTOOLS_INDEX_ANNOTATED(ch_annotate)

    ch_stats_input = ch_annotate
    .combine(ch_annotate_index, by: 0)

    BCFTOOLS_STATS(ch_stats_input)

    }

}