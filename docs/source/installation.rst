.. _installation:

Installation
=======================

Dependencies
-----------------------

* `libhdf5 <https://support.hdfgroup.org/HDF5/release/obtain5.html>`_ is automatically downloaded and compiled when running make but this can be disabled with: HDF5=nofetch make. The nanopolish binary will link libhdf5.a statically.
* `eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ is also automatically downloaded and included when compiling with make.
* `biopython <http://biopython.org/>`_ is required to run the helpers in scripts/.
* `htslib <https://github.com/samtools/htslib>`_ is included as a submodule and compiled automatically.
* A compiler that supports C++11 is needed to build the sources. Development of the code is performed using gcc-4.8.

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
