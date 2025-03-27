#!/bin/bash

SAC_DISPLAY_COPYRIGHT=0

out="compare_pdf"
rm -rf $out
mkdir -p $out 


for synpath in syn*; do 
for cout in $synpath/cout*; do
fname=$(basename $cout)
echo $cout $(dirname $cout)/p${fname:1:}
sac << EOF
r $cout $(dirname $cout)/p${fname:1}
qdp off
color red inc list blue red 
p1
saveimg $out/compare_${synpath}_${fname}.pdf
q
EOF

done
done


pdftk $out/* cat output out.pdf