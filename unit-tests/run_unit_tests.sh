#!/bin/bash

dir="$(dirname $0)"

echo "Run all unit tests ..."

# 7 April 2021
# https://stackoverflow.com/questions/2437452/how-to-get-the-list-of-files-in-a-directory-in-a-shell-script
for file in "$dir"/*; do
    if [[ "$file" != "$0" ]]; then
        $file
    fi
done

echo "Done."
