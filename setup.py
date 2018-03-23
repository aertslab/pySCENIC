# -*- coding: utf-8 -*-
"""A setuptools based setup module.

For more information, please see:
- https://pypi.python.org/pypi/setuptools
- https://pythonhosted.org/setuptools
- https://python-packaging-user-guide.readthedocs.io/en/latest/distributing/
- https://packaging.python.org/en/latest/distributing.html
- https://github.com/pypa/sampleproject

"""
import glob
import io
import os
import setuptools
import versioneer
from os.path import join


def read(*names, **kwargs):
    r"""Return the contents of a file.

    Default encoding is UTF-8, unless specified otherwise.

    Args:
        - names (list, required): list of strings, parts of the path.
          the path might be relative to the current file.

    Keyword Args:
        **kwargs: Arbitrary keyword arguments.

    Returns:
      String containing the content of the file.

    Examples:
        >>> read('docs/readme.rst')
            u'Overview\n--------\n...'

        >>> read('docs', 'readme.rst')
            u'Overview\n--------\n...'

    """
    fn = os.path.join(os.path.dirname(__file__), *names)
    with io.open(fn, encoding=kwargs.get('encoding', 'utf8')) as fd:
        return fd.read()


def read_requirements(fname):
    with open(fname, 'r', encoding='utf-8') as file:
        return [line.rstrip() for line in file]


setuptools.setup(
    name='pyscenic',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Python implementation of the SCENIC pipeline for transcription factor inference from single-cell transcriptomics experiments.",
    long_description=read('README.rst'),
    classifiers=[
        # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    python_requires=">=3.5",
    keywords='single-cell transcriptomics gene-regulatory-network transcription-factors',
    author="Bram Van de Sande",
    url='https://github.com/aertslab/pySCENIC',
    license='GPL-3.0+',
    packages=setuptools.find_packages(where='src'),
    package_dir={'': 'src'},
    py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob.glob('src/*.py')],
    include_package_data=True,
    install_requires=read_requirements('requirements.txt'),
    entry_points = {
        'console_scripts': ['pyscenic = pyscenic.scenic:scenic',
                            'db2feather = pyscenic.db2feather:main'],
    }
)