#!/bin/bash 

dist=3

depsrc=3
deprcv=0

nt=500
dt=0.01

# ----------------------- model with soil ---------------------
modname="killari"
out="GRN"

# compute Green's Functions
grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# generate synthetic seismograms
grt.syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -A45 -S1e24 -M0/45/90 -D0/tfunc -Osyn


# ----------------------- model without soil ---------------------
modname="killari_nosoil"
out="GRN_nosoil"

# compute Green's Functions
grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# generate synthetic seismograms
grt.syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -A45 -S1e24 -M0/45/90 -D0/tfunc -Osyn_nosoil

