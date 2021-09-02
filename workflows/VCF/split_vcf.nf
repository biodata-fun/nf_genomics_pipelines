#!/usr/bin/env nextflow

/* 
 * Script for splitting a VCF containing multiple chromosomes into 1 VCF per chromosome.
 * It will only consider the biallelic sites
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

nextflow.enable.dsl=2

