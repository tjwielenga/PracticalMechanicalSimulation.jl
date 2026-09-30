PANDOC ?= pandoc
SOURCE := sparseFullyConsistentMethod.md
OUTPUT := sparse-fully-consistent-modeling-method.docx
EQUATION_FILTER := equation-numbering.lua

.PHONY: docx

docx:
	$(PANDOC) $(SOURCE) --from=markdown-implicit_figures \
		--lua-filter=$(EQUATION_FILTER) --citeproc --resource-path=. -o $(OUTPUT)
