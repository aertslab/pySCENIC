#!/usr/bin/env bash

# Install pyenv:
#
# https://github.com/pyenv/pyenv#automatic-installer
#
# Install pyenv in ~/.pyenv:
#
#     git clone https://github.com/pyenv/pyenv.git ~/.pyenv
#
# Optionally, try to compile a dynamic Bash extension to speed up Pyenv.
#
#     cd ~/.pyenv && src/configure && make -C src
#
# Get list of possible python interpreters to install:
#    ~/.pyenv/bin/pyenv install -l
#
# Install latest versions of python interpreters
#
#    ~/.pyenv/bin/pyenv install 3.7.13
#    ~/.pyenv/bin/pyenv install 3.8.13
#    ~/.pyenv/bin/pyenv install 3.9.13
#    ~/.pyenv/bin/pyenv install 3.10.6
#
# Enable those versions for pySCENIC (same dir as setup.py):
#    ~/.pyenv/bin/pyenv local 3.7.13 3.8.13 3.9.13 3.10.6

# Add pyenv to PATH:
export PYENV_ROOT="${HOME}/.pyenv"
command -v pyenv >/dev/null || export PATH="${PYENV_ROOT}/bin:${PATH}"

eval "$(pyenv init -)"
tox
