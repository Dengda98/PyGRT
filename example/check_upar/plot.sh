#!/bin/bash

SAC_DISPLAY_COPYRIGHT=0

out="compare_pdf"
rm -rf $out
mkdir -p $out

diff=0.001
depsrc=2
deprcv=3.2
deprcv2=`echo $deprcv+$diff | bc | awk '{printf("%.3f", $1)}'`
dist=10
dist2=`echo $dist+$diff | bc | awk '{printf("%.3f", $1)}'`

GRN="GRN"

# --------------------- ui_z --------------------------------------
for ch in $(ls $GRN/milrow_${depsrc}_${deprcv2}_${dist}); do
echo $ch
sac <<EOF
r $GRN/milrow_${depsrc}_${deprcv2}_${dist}/$ch
subf $GRN/milrow_${depsrc}_${deprcv}_${dist}/$ch
div -$diff
w uiz_diff.sac
r uiz_diff.sac
r more $GRN/milrow_${depsrc}_${deprcv}_${dist}/z${ch:0:3}.sac
int
qdp off
color red inc list blue red 
p2
saveimg $out/compare_$ch.pdf
q
EOF

done

rm uiz_diff.sac

pdftk $out/* cat output compare_uiz.pdf

rm -rf $out
mkdir -p $out

# --------------------- ui_z --------------------------------------
for ch in $(ls $GRN/milrow_${depsrc}_${deprcv}_${dist2}); do
echo $ch
sac <<EOF
r $GRN/milrow_${depsrc}_${deprcv}_${dist2}/$ch
subf $GRN/milrow_${depsrc}_${deprcv}_${dist}/$ch
div $diff
w uir_diff.sac
r uir_diff.sac
r more $GRN/milrow_${depsrc}_${deprcv}_${dist}/r${ch:0:3}.sac
int
qdp off
color red inc list blue red 
p2
saveimg $out/compare_$ch.pdf
q
EOF

done

rm uir_diff.sac

pdftk $out/* cat output compare_uir.pdf

rm -rf $out