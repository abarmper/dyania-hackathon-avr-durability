#!/usr/bin/env bash
# Build report/main.pdf. Needs TeX Live with biber (present at /usr/bin).
cd "$(dirname "$0")" && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex >/dev/null && echo "built main.pdf"
