from setuptools import find_packages, setup

setup(name='kython',
packages=find_packages(include=['kython']),
version='0.1',
description='Python Library to compute and plot phylogeny Kmer signatures',
author='Adrien Leroy',
license='MIT',
install_requires=[],
setup_requires=['pytest-runner'],
tests_require=['pytest'],
test_suite='tests')
