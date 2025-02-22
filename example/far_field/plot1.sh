#!/bin/bash

sac <<EOF
r GRN/milrow_2_0_1800/EXZ.sac
r more GRN_filon/milrow_2_0_1800/EXZ.sac
int
qdp off
color red inc list blue red 
width 1
# p2
p1
save FIM_compare1.pdf
q
EOF