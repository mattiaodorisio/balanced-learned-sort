#!/bin/bash

cd data

# Remove files without extension that are not directory
ls -p -I '*.*' | grep -v / | xargs rm

# Remove .txt files in subdirectories
rm */*.txt

cd ..
