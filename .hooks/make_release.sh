#!/bin/bash

# 14 May 2023
# https://stackoverflow.com/a/16496491

usage() {
    echo "Usage: $0 -v x.x.x " 1>&2;
    exit 1
}

while getopts ":v:" o; do
    case "${o}" in
        v)
            new_version=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${new_version}" ]; then
    usage
fi

fname="../configure.ac"

# get current version from configure.ac
current_version=$(sed -nre 's/AC_INIT\(\[epic\], \[([0-9]+\.[0-9]+\.[0-9]+)\], (.*)\)/\1/p' ${fname})

echo "Update version string from ${current_version} to ${new_version}"

sed -i "s:${current_version}:${new_version}:g" ${fname}

git add ${fname}
git commit -m "Update version string from ${current_version} to ${new_version}"
git push
