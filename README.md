# Nanopolish

A nanopore consensus algorithm using a signal-level hidden Markov model.

## Dependencies

The program requires [libhdf5](http://www.hdfgroup.org/HDF5/release/obtain5.html) and a compiler that supports C++11. Development of the code is performed using [gcc-4.8](https://gcc.gnu.org/gcc-4.8/). libhdf5 can be automatically installed by the Makefile if you do not have it already (see below).

## Installation instructions

You will need to run ```git clone --recursive https://github.com/jts/nanopolish.git``` to get the source code and submodules. You can then compile nanopolish by running:

```
make
```

If you do not already have libhdf5, you can tell nanopolish to download and build it locally in the source tree:

```
make libhdf5.install nanopolish
```

## Brief usage instructions

The pipeline is still a prototype so it is fragile at the moment. It will be revised for general use after we submit the paper.

The reads that are input into the HMM must be output as a ```.fa``` file  by ```poretools```. This is important as ```poretools``` writes the path to the original ```.fast5``` file (containing the signal data) in the fasta header. These paths must be correct or nanopolish cannot find the events for each read. Let's say you have exported your reads to ```reads.fa``` and you want to polish ```draft.fa```. You should run:

```
make -f consensus.make READS=reads.fa ASSEMBLY=draft.fa
```

This will map the reads to the assembly with ```bwa mem -x ont2d``` and export a file mapping read names to fast5 files.

You can then run ```nanopolish consensus```. It is recommended that you run this in parallel.

```
python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 nanopolish consensus -o nanopolish.{1}.fa -w {1} --r reads.pp.fa -b reads.pp.sorted.bam -g draft.fa -t 4
```

This command will run the consensus algorithm on eight 100kbp segments of the genome at a time, using 4 threads each. Change the ```-P``` and ```--threads``` options as appropriate for the machines you have available.

After all polishing jobs are complete, you can merge the individual segments together into the final assembly:

```
python nanopolish_merge.py draft.fa nanopolish.*.fa > polished.fa
```

## Known Issues

If you have extremely high depth (for example you sequenced a virus) then ```poa``` and the hmm will take a very long time to run. I suggest downsampling to reasonable coverage before trying to call the consensus sequence.

## Credits and Thanks

The fast table-driven logsum implementation was provided by Sean Eddy as public domain code. This code was originally part of [hmmer3](http://hmmer.janelia.org/).
