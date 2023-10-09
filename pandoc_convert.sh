#!/bin/bash
pandoc -t markdown-citations-simple_tables -crossref -s supplements.tex -o supplements.md --bibliography bibliography.bib --citeproc --wrap=preserve
pandoc -t markdown-citations-simple_tables -crossref -s main.tex -o main.md --bibliography bibliography.bib --citeproc --wrap=preserve