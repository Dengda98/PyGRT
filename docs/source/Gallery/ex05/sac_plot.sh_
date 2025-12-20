#!/bin/bash

SAC_DISPLAY_COPYRIGHT=0

out="compare_pdf"
rm -rf $out
mkdir -p $out 


for synpath in syn_ex syn_sf syn_dc syn_mt; do 
for ch in Z R T; do
cname=c_$synpath/$ch.sac
pname=p_$synpath/$ch.sac
echo $cname $pname
sac << EOF
r $cname $pname
qdp off
color red inc list blue red 
p1
saveimg $out/compare_${synpath}_${ch}.pdf
q
EOF

done
done


pdftk $out/* cat output out.pdf