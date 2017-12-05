.. _installation:

Installation
=======================

Dependencies
-----------------------

A compiler that supports C++11 is needed to build nanopolish. Development of the code is performed using [gcc-4.8](https://gcc.gnu.org/gcc-4.8/).

By default, nanopolish will download and compile all of its required dependencies. Some users however may want to use system-wide versions of the libraries. To turn off the automatic installation of dependencies set `HDF5=noinstall`, `EIGEN=noinstall` or `HTS=noinstall` parameters when running `make` as appropriate. The current versions and compile options for the dependencies are:

* `libhdf5-1.8.14 <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_ compiled with multi-threading support ``--enable-threadsafe``
* `eigen-3.2.5 <http://eigen.tuxfamily.org/>`_ 
* `htslib-1.4 <http://github.com/samtools/htslib>`_

Additionally the helper `scripts` require `biopython <http://biopython.org/>`_ and `pysam <http://pysam.readthedocs.io/en/latest/installation.html>`_.

Installing the latest code from github (recommended)
------------------------------------------------------
You can download and compile the latest code from github as follows ::

    git clone --recursive https://github.com/jts/nanopolish.git
    cd nanopolish
    make

Installing a particular release
------------------------------------------------------
When major features have been added or bugs fixed, we will tag and release a new version of nanopolish. If you wish to use a particular version, you can checkout the tagged version before compiling ::

    git clone --recursive https://github.com/jts/nanopolish.git
    cd nanopolish
    git checkout v0.7.1
    make

To run using docker
-------------------

First build the image from the dockerfile: ::

    docker build .

Note the uuid given upon successful build. Then you can run nanopolish from the image: ::

    docker run -v /path/to/local/data/data/:/data/ -it :image_id  ./nanopolish eventalign -r
 /data/reads.fa -b /data/alignments.sorted.bam -g /data/ref.fa
