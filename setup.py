from setuptools import find_packages, setup

setup(name='kython',
packages=find_packages(include=['kython']),
version='1.0.3',
description='Python Library to compute and plot phylogeny Kmer signatures',
author='Adrien Leroy, Alex Lence, Simon Chardin',
author_email='leroy.adrien.76@gmail.com',
long_description="University project to understand how to build python libraries while trying to compute and plot Kmer signature from whole genome on the ncbi website",
url="http://github.com/LeroyAdrien/Kython",
license='MIT',
install_requires=['biopython','numpy','matplotlib','pandas','ncbi-genome-download','scipy'],
setup_requires=['pytest-runner'],
tests_require=['pytest'],
test_suite='tests')
