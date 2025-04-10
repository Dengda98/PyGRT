#!/bin/bash 


mod="halfspace2"
stgrt -M$mod -D20/16 -X3/3/1 -Y-10/10/41 -e > grn1

stgrt.syn -M0/70/0 -Su1e6 -e -N < grn1 > syn_1_stkslip
stgrt.stress < syn_1_stkslip > stress_1_stkslip

stgrt.syn -M0/70/90 -Su1e6 -e -N < grn1 > syn_1_dipslip
stgrt.stress < syn_1_dipslip > stress_1_dipslip



stgrt -M$mod -D10/14 -X3/3/1 -Y-10/10/41 -e > grn2

stgrt.syn -M0/70/0 -Su1e6 -e -N < grn2 > syn_2_stkslip
stgrt.stress < syn_2_stkslip > stress_2_stkslip

stgrt.syn -M0/70/90 -Su1e6 -e -N < grn2 > syn_2_dipslip
stgrt.stress < syn_2_dipslip > stress_2_dipslip


mod="mod"
stgrt -M$mod -D5/0 -X3/3/1 -Y-10/10/41 -e > grn3

stgrt.syn -M0/70/0 -Su1e6 -e -N < grn3 > syn_3_stkslip
stgrt.stress < syn_3_stkslip > stress_3_stkslip

stgrt.syn -M0/70/90 -Su1e6 -e -N < grn3 > syn_3_dipslip
stgrt.stress < syn_3_dipslip > stress_3_dipslip
