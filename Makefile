# This is a comment.
# You must indent using <TAB>s, not spaces.

# If you're using GNU make and you need help debugging a makefile then there
# is a single line your should add. And it's so useful that you should add it
# to every makefile you create.
# cfr. http://blog.jgc.org/2015/04/the-one-line-you-should-add-to-every.html
print-%: ; @echo $*=$($*)


# notdir: Extracts all but the directory-part of each filename.
# cfr https://www.gnu.org/software/make/manual/html_node/File-Name-Functions.html
PACKAGENAME = $(notdir $(CURDIR))

# For ENVNAME, we assume the package-code is located in a 'src'-subdirectory of the virtualenv.
ENVNAME = $(notdir $(realpath $(CURDIR)/../.. ))


clean:
	- rm -rf dist/
	- rm -fr docs/_build

flake8:
	- flake8 src/$(PACKAGENAME)

py:
	- py.test

py-pdb:
	- py.test --pdb

py-cov:
	- py.test --cov=$(PACKAGENAME) --cov-report html -v tests/

py-pdb-cov:
	- py.test --cov=$(PACKAGENAME) --cov-report html -v tests/  --pdb

py-pdb-cov-lf:
	- py.test --cov=$(PACKAGENAME) --cov-report html -v tests/  --pdb --lf
