# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

# export TEXINPUTS=$PWD"/source/_latex/:"

#macro
REMOVEAPPLEJUNK = find source/ -iname "._*" -exec rm {} \;

.PHONY: help clean html web pickle htmlhelp latex changes linkcheck

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  upload    to upload to Sourceforge"
	@echo "  sitemap   update the google sitemap"
	@echo "  all       clean, html, latex, sitemap, upload"


all: html latex sitemap upload

clean:
	-rm -rf build/*
	$(REMOVEAPPLEJUNK)

html:
	$(REMOVEAPPLEJUNK)
	mkdir -p build/html build/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) build/html
	@echo
	@echo "Build finished. The HTML pages are in build/html."

pickle:
	$(REMOVEAPPLEJUNK)
	mkdir -p build/pickle build/doctrees
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) build/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

web: pickle

json:
	$(REMOVEAPPLEJUNK)
	mkdir -p build/json build/doctrees
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) build/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(REMOVEAPPLEJUNK)
	mkdir -p build/htmlhelp build/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in build/htmlhelp."

latex:
	$(REMOVEAPPLEJUNK)
	mkdir -p build/latex build/doctrees
	cp source/_latex/*.* build/latex/
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex
	@echo
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."
	cd build/latex; make all-pdf

changes:
	$(REMOVEAPPLEJUNK)
	mkdir -p build/changes build/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) build/changes
	@echo
	@echo "The overview file is in build/changes."

linkcheck:
	$(REMOVEAPPLEJUNK)
	mkdir -p build/linkcheck build/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/linkcheck/output.txt."

sitemap:
	cd sitemap_gen; python sitemap_gen.py --config=config.xml
	
	
upload:
	@echo "Uploading to SourceForge.net..."
	$(REMOVEAPPLEJUNK)
	rsync -rtza --rsh=ssh build/html/* rhw27,rmg@web.sourceforge.net:htdocs/3.2/
	rsync -rtza --rsh=ssh build/latex/*.pdf rhw27,rmg@web.sourceforge.net:htdocs/3.2/
	#scp -r build/html/* athena.dialup.mit.edu:www/RMG/
	#scp build/latex/*.pdf athena.dialup.mit.edu:www/RMG/