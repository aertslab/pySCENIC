#!/usr/bin/env bash

# Installation for MacOS X
# 1. brew update
# 2. brew install pyenv
# Check latest versions of python interpreters on
# 3. pyenv install 3.5.5 && pyenv install 3.6.4
# In same directory as setup.py:
# 4. pyenv local 3.5.5 3.6.4

eval "$(pyenv init -)"
tox
