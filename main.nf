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

    if ( params.udg ) {
        ch_input = BAM2SNPAD_UDG(ch_bam, ch_bai, ch_ref_fasta)
    } else ( !params.udg ) {
        ch_input = BAM2SNPAD_NOT_UDG(ch_bam, ch_bai, ch_ref_fasta)
    }

    ch_map_bed = Channel.from(params.map_bed)

    MAPABILITY(ch_input, ch_map_bed)
    ch_snpad_input = MAPABILITY.out

    SNPAD(ch_snpad_input)
    ch_snpad_params = SNPAD.out.inputs

    CALL(ch_snpad_params, MAPABILITY.out)
    ch_split_vcfs = CALL.out

}