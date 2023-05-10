#!/bin/bash

dir="$(dirname $0)"
nprocs=4

echo "Run all unit tests ..."

# 7 April 2021
# https://stackoverflow.com/questions/2437452/how-to-get-the-list-of-files-in-a-directory-in-a-shell-script
for file in "$dir"/*; do
    if [[ "$file" != "$0" ]]; then
        # 14 April 2022
        # https://stackoverflow.com/questions/229551/how-to-check-if-a-string-contains-a-substring-in-bash
        if [[ "$file" == *"mpi"* ]]; then
            mpirun -np $nprocs $file
        else
            $file
        fi
    fi
done

echo "Done."
