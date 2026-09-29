PANDOC ?= pandoc
SOURCE := sparseFullyConsistentMethod.md
OUTPUT := sparse-fully-consistent-modeling-method.docx

.PHONY: docx

docx:
	$(PANDOC) $(SOURCE) --from=markdown-implicit_figures --citeproc --resource-path=. -o $(OUTPUT)
