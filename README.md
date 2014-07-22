Grow-Cut Image Segmentation
==========================

The "growcut" package implements cellular automata based segmentation algorithm.

References:
-----------

TL Kline, M Zamir, EL Ritman, 
"Accuracy of microvascular measurements obtained from micro-CT images,"
Annals of biomedical engineering 38 (9), 2851-2864, 2010


Maintainers
-----------

    - Timothy Lee Kline

Contributors
------------
    - 

Dependencies
------------

The required dependencies to build the software are:

  - python
  - numpy
  - scipy
  - skfmm
  - nibabel

Install
-------

This packages uses distutils, which is the default way of installing python modules. To install in your home directory, use:

    python setup.py install --home

To install for all users on Unix/Linux:

    python setup.py build
    sudo python setup.py install

Development
-----------

Follow: Fork + Pull Model:

    https://help.github.com/articles/using-pull-requests