#!/bin/bash
echo "Compiling FAUST file"
faust -mdoc -mdlang fr $1
inkscape --export-type pdf "$(basename -- "$1" .dsp)-mdoc"/svg/svg-01/process.svg
pdflatex -output-directory "$(basename -- "$1" .dsp)-mdoc"/pdf "$(basename -- "$1" .dsp)-mdoc"/tex/"$(basename -- "$1" .dsp).tex"
pdflatex -output-directory "$(basename -- "$1" .dsp)-mdoc"/pdf "$(basename -- "$1" .dsp)-mdoc"/tex/"$(basename -- "$1" .dsp).tex"