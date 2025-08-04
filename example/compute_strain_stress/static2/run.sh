#!/bin/bash 


mod="halfspace2"
grt static greenfn -M$mod -D20/16 -X3/3/1 -Y-10/10/0.5 -e > grn1

grt static syn -M0/70/0 -Su1e6 -e -N < grn1 > syn_1_stkslip
grt static stress < syn_1_stkslip > stress_1_stkslip

grt static syn -M0/70/90 -Su1e6 -e -N < grn1 > syn_1_dipslip
grt static stress < syn_1_dipslip > stress_1_dipslip



grt static greenfn -M$mod -D10/14 -X3/3/1 -Y-10/10/0.5 -e > grn2

grt static syn -M0/70/0 -Su1e6 -e -N < grn2 > syn_2_stkslip
grt static stress < syn_2_stkslip > stress_2_stkslip

grt static syn -M0/70/90 -Su1e6 -e -N < grn2 > syn_2_dipslip
grt static stress < syn_2_dipslip > stress_2_dipslip


mod="mod"
grt static greenfn -M$mod -D5/0 -X3/3/1 -Y-10/10/0.5 -e > grn3

grt static syn -M0/70/0 -Su1e6 -e -N < grn3 > syn_3_stkslip
grt static stress < syn_3_stkslip > stress_3_stkslip

grt static syn -M0/70/90 -Su1e6 -e -N < grn3 > syn_3_dipslip
grt static stress < syn_3_dipslip > stress_3_dipslip
