# -*- coding: utf-8 -*-

from fabric.api import run, env, local, task, put
from fabric.context_managers import lcd, cd


# Use servernames from ~/.ssh/config file.
env.use_ssh_config = True
# List of remote serves.
env.hosts = ['hpc2-big7']


@task
def clean():
    with lcd('..'):
        local('rm -rf dist/*')
    with cd('~/data/downloads/'):
        run('rm -rf *')


@task
def deploy(version):
    clean()

    with lcd('..'):
        local('git tag -a {0} -m \"version {0}\"'.format(version))
        local('git push --tags')

        local('python setup.py sdist')
        #local('python setup.py bdist_wheel')
        put('./dist/pyscenic*', '~/data/downloads')

    with cd('~/data/downloads/'):
        run('. activate pyscenic')
        run('pip install --upgrade pyscenic-*.tar.gz')


@task
def install():
    # Install miniconda
    # Create conda environment
    # Install python package pyscenic in environment
    pass


