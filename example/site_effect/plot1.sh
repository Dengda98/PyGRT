#!/bin/bash

sac <<EOF
r syn/*
qdp off
p1
save syn1.pdf

r syn_nosoil/*
qdp off
p1
save syn1_nosoil.pdf
q
EOF