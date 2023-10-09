#!/bin/bash
quarto pandoc -t markdown-citations-simple_tables -crossref -s supplements.tex -o supplements.qmd --bibliography bibliography.bib --citeproc --wrap=preserve
quarto pandoc -t markdown-citations-simple_tables -crossref -s main.tex -o main.qmd --bibliography bibliography.bib --citeproc --wrap=preserve