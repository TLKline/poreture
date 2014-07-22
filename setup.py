""" poreture package configuration """

from setuptools import setup

setup(
    name='poreture',
    version='0.1',
    description='Poreture - Vasculature and Pore Network Analysis Toolkit.',
    author='Timothy Lee Kline',
    author_email='kline.timothy@mayo.edu',
    packages=['poreture'],
    install_requires=[
        'nibabel',
	    'numpy',
        'scipy'])