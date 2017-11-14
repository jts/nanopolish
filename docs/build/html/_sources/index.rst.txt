.. Nanopolish documentation master file, created by
   sphinx-quickstart on Tue Nov 14 14:59:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Nanopolish
======================================
`Nanopolish <https://github.com/jts/nanopolish>`_ is a software package for signal-level analysis of Oxford Nanopore sequencing data. Nanopolish can calculate an improved consensus sequence for a draft genome assembly, detect base modifications, call SNPs and indels with respect to a reference genome and more (see Nanopolish modules, below).


.. toctree::
   :hidden:

   installation
   quickstart
   manual
   faq

Publication
====================================
Loman, Nicholas J., Joshua Quick, and Jared T. Simpson. "A complete bacterial genome assembled de novo using only nanopore sequencing data." Nature methods 12.8 (2015): 733-735.

Credits and Thanks
========================
The fast table-driven logsum implementation was provided by Sean Eddy as public domain code. This code was originally part of `hmmer3 <http://hmmer.org/>`_ . Nanopolish also includes code from Oxford Nanopore's `scrappie <https://github.com/nanoporetech/scrappie>`_ basecaller. This code is licensed under the MPL.
