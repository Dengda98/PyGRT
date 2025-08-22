#!/bin/bash


# No Upsampling
grt greenfn -Mmilrow -N200/0.5 -D1/0 -R100 -OGRN

# Upsampling 5 times
grt greenfn -Mmilrow -N200/0.5+n5 -D1/0 -R100 -OGRN_up