#!/bin/env bash

version=$(sed -nre 's/AC_INIT\(\[epic\], \[([0-9]+\.[0-9]+\.[0-9]+)\], (.*)\)/\1/p' ../configure.ac)

cd manual

sed -i "s:@VERSION@:$version:g" intro.adoc

asciidoctor intro.adoc --destination-dir=../html
asciidoctor-pdf intro.adoc --destination-dir=../pdf

sed -i "s:$version:@VERSION@:g" intro.adoc
