#!/bin/bash
quarto pandoc -t markdown -s supplements.tex -o supplements.qmd --bibliography bibliography.bib --wrap=preserve
quarto pandoc -t markdown -s main.tex -o main.qmd --bibliography bibliography.bib --wrap=preserve