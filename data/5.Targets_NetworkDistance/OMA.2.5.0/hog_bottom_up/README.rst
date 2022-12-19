GETHOGs: A Tool to Infer Hierarchical Orthologous Groups (HOGs)
===============================================================


Motivation
----------
GETHOGs is a library to infer hierarchical orthologous groups from pairwise evolutionary relationship and its related reference species tree.
Infered HOGs are provided as orthoXML ( see http://orthoxml.org.). GETHOGs can be used as a python library or using the CLI script in the bin folder.

For extended documentation, please check the related paper "Clément-Marie Train, Natasha M. Glover, Gaston H. Gonnet, Adrian M. Altenhoff, Christophe Dessimoz; Orthologous Matrix (OMA) algorithm 2.0: more robust to asymmetric evolutionary rates and more scalable hierarchical orthologous group inference, Bioinformatics, Volume 33, Issue 14, 15 July 2017, Pages i75–i82, https://doi.org/10.1093/bioinformatics/btx229"
or look at the docs folder that contain information
on common use cases and API documentation of the library.


Installation
------------
GETHOGs is written in python3.

Consider using pip to install the package directly from a checked out git repo

.. code-block:: sh

   python -m pip install --upgrade pip
   pip install gethogs

