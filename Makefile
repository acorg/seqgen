.PHONY: pytest, check, tcheck, pycodestyle, pyflakes, flake8, lint, wc, clean, clobber, upload

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)
PYDIRS := bin seqgen test

pytest:
	pytest

check:
	python -m discover -v

tcheck:
	trial --rterrors test/test*py

pycodestyle:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pycodestyle

pyflakes:
	find .  -path './.tox' -prune -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pyflakes

flake8:
	find $(PYDIRS) -name '*.py' -print0 | $(XARGS) -0 flake8 --ignore E402,W504

lint: pycodestyle pyflakes

wc:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	find . -name '__pycache__' -type d -print0 | xargs -0 rmdir
	find . -name '_trial_temp' -type d -print0 | xargs -0 rm -r
	python setup.py clean

clobber: clean
	rm -fr .tox seqgen.egg-info dist

# The upload target requires that you have access rights to PYPI. You'll also
# need twine installed (on OS X with brew, run 'brew install twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/seqgen-$$(env PYTHONPATH=.:$$PYTHONPATH bin/seq-gen-version.py).tar.gz
