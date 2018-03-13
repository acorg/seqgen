#!/usr/bin/env python

from setuptools import setup


# Modified from http://stackoverflow.com/questions/2058802/
# how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package
def version():
    import os
    import re

    init = os.path.join('seqgen', '__init__.py')
    with open(init) as fp:
        initData = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]+)['\"]",
                      initData, re.M)
    if match:
        return match.group(1)
    else:
        raise RuntimeError('Unable to find version string in %r.' % init)


setup(name='seqgen',
      version=version(),
      packages=['seqgen'],
      url='https://github.com/acorg/seqgen',
      download_url='https://github.com/acorg/seqgen',
      author='Terry Jones',
      author_email='tcj25@cam.ac.uk',
      keywords=['genetic sequences'],
      classifiers=[
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      description=('Command line script and Python class for generating '
                   'genetic sequence data in FASTA format.'),
      long_description=(
          'Please see https://github.com/acorg/seqgen for details.'),
      license='MIT',
      scripts=['bin/seq-gen.py', 'bin/seq-gen-version.py'],
      install_requires=['dark-matter>=1.1.28'])
