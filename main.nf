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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BAM2SNPAD_NOT_UDG } from './modules/local/bam2snpad/not_udg/main'
include { BAM2SNPAD_UDG     } from './modules/local/bam2snpad/udg/main'
include { MAPABILITY        } from './modules/local/mapability/main'
include { SNPAD             } from './modules/local/snpAD/main'
include { CALL              } from './modules/local/call/main'
include { CONCAT            } from './modules/local/concat/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


workflow {

    ch_bam = Channel.fromFilePairs(params.inputbam, size: -1)
        .map {
            meta, bam ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            [ fmeta, bam ]
        }

    ch_bai = Channel.fromPath(params.bai)
    ch_ref_fasta = Channel.from(params.ref_fasta)
    //ch_i = Channel.of(1..22, 'X', 'Y')
    ch_i = Channel.of(1..3)

    ch_input_bam2snpad = ch_bam
    .combine(ch_bai)
    .combine(ch_ref_fasta)
    .combine(ch_i)
    .multiMap{
            meta, bam, bai, fasta, chrom ->
            bam:      [ meta, bam ]
            bai:      [ bai       ]
            fasta:    [ fasta     ]
            chrom:    [ chrom     ]
        }
    
    if ( params.udg ) {
        ch_input = BAM2SNPAD_UDG(ch_input_bam2snpad.bam, ch_input_bam2snpad.bai, ch_input_bam2snpad.fasta, ch_input_bam2snpad.chrom.flatten())
    } else ( !params.udg ) {
        ch_input = BAM2SNPAD_NOT_UDG(ch_input_bam2snpad.bam, ch_input_bam2snpad.bai, ch_input_bam2snpad.fasta, ch_input_bam2snpad.chrom.flatten())
    }


    ch_map_bed = Channel.from(params.map_bed)

    //ch_input_mappability =
    ch_input
    .combine(ch_map_bed)
    .multiMap{
            meta, input, chrom, map ->
            input:    [ meta, input ]
            chrom:    [ chrom       ]
            map:      [ map         ]
        }

    // ch_input_mappability.input.view()

    // MAPABILITY(ch_input, ch_map_bed)
    // ch_snpad_input = MAPABILITY.out

    // SNPAD(ch_snpad_input)
    // ch_snpad_params = SNPAD.out

    // ch_ref_fasta_fai = Channel.from(params.ref_fasta_fai)
    // CALL(ch_snpad_params, MAPABILITY.out, ch_ref_fasta_fai)
    // ch_split_vcfs = CALL.out

    // ch_rsID = Channel.from(params.rsID)
    // CONCAT(ch_split_vcfs, params.rsID)

}