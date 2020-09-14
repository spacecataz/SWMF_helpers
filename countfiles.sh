#!/bin/bash

for f in $@; do 
    echo $f
    ls -aR $f | wc -l
done
