#!/bin/bash 

sac<<EOF 
r GRN/*/VFZ.sac
sss
tw 0 15 
prs size 0.3 l off
saveimg view.pdf
q
EOF