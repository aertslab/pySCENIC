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

long_description = (
        read('docs/readme.rst') +
        '\n\n' +
        read('docs/changes.rst') +
        '\n\n' +
        read('docs/contributors.rst')
)

setuptools.setup(
    name='pyscenic',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Python implementation of the SCENIC pipeline.",
    long_description=long_description,
    classifiers=[
        # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    keywords='single-cell rna-seq',
    author="Bram Van de Sande",
    author_email="vandesande.bram@gmail.com",
    url='http://pyscenic.readthedocs.io/en/latest/',
    license='GPL-3.0+',
    packages=setuptools.find_packages('pyscenic'),
    package_dir={'': 'pyscenic'},
    py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob.glob('pyscenic/*.py')],

    include_package_data=True,
    zip_safe=False,
    #test_suite='nose.collector',

    install_requires=[
        'dask',
    ],

    setup_requires=[
        'pytest-runner'
    ],

    tests_require=[
        'coverage',
        'flake8',
        'flake8-blind-except',
        'flake8-coding',
        'flake8-docstrings',
        'flake8-exact-pin',
        'flake8-import-order',
        'flake8-pep3101',
        'flake8-print',
        'flake8-respect-noqa',
        #'flask-testing',
        'mock',
        'pook',
        'pytest',
        'pytest-capturelog',
        'pytest-cov',
        'pytest-flake8',
        'pytest-mccabe',
        'pytest-mock',
    ],

    scripts=[
        'scripts/flaskapp_ucb_tas.py',
        'scripts/uwsgitop.sh',
    ],

    entry_points={
        'console_scripts': [
            'ucb-tas = ucb_tas.cli.run:run',
            'ucb-tas-celery = ucb_tas.cli.run:run_celery',
            'ucb-tas-check-system = ucb_tas.cli.cli_non_app:check_system',
        ],

    },
)