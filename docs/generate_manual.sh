#!/bin/env bash

cd manual
asciidoctor intro.adoc --destination-dir=../html
asciidoctor-pdf intro.adoc --destination-dir=../pdf
