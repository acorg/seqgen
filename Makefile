.PHONY: check, pycodestyle, pyflakes, wc, clean, clobber, upload

check:
	pytest

pycodestyle:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pycodestyle

pyflakes:
	find .  -path './.tox' -prune -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pyflakes

wc:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	find . -name '__pycache__' -type d -print0 | xargs -0 rmdir
	find . -name '_trial_temp' -type d -print0 | xargs -0 rm -r
	find . -name '.mypy_cache' -type d -print0 | xargs -0 rm -r
	find . -name '.pytest_cache' -type d -print0 | xargs -0 rm -r
	rm -fr seqgen.egg-info
	python setup.py clean

clobber: clean
	rm -fr .tox seqgen.egg-info dist

# The upload target requires that you have access rights to PYPI. You'll also
# need twine installed (on OS X with brew, run 'brew install twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/seqgen-$$(env PYTHONPATH=.:$$PYTHONPATH bin/seq-gen-version.py).tar.gz
