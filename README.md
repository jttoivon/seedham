Installing and pre-requisities
==============================

As the only dependency, seedham requires Boost library, which is normally installed on Linux machines. Relatively recent version
should be used. At least version 1.49 is known to work.
Running make in the directory of the distribution should compile 'seedham'.
If wanted, you can install by running the command

	sudo make install

or if you want to install to a non-standard location use, for example

	make prefix=$HOME/usr install

which installs the executables (multinomial, seedham, seedham+) to $HOME/usr/bin.
Running command 'seedham' or 'seedham+' should give brief instructions on the command line parameters.
The distribution also includes, for internal use, an implementation of suffix array by Juha Kärkkäinen in directory CPM03.

Running
=======

The input file to seedham should consist of sequences separated by new lines. If you give a fasta file as
parameter, it ignores the header lines (ones that begin with '>'). If a fasta sequence is split on
several lines, then seedham considers these as separate sequences. So, using a fasta file at the moment
is not recommended. Also, currently sequences containing non-base characters, such as 'N', are ignored.

Examples of running seedham:

	./seedham 8 data/TFAP2A-head-1000.seq            # By default, Hamming radius 1 is used. Use best seed of length 8

	./seedham -r 3 GGGCA data/TFAP2A-head-1000.seq   # Explicitly set Hamming radius to be 3. Use seed GGGCA

The program seedham+ has the same command line parameters as seedham.

The data file data/TFAP2A-head-1000.seq contains the first 1000 reads from ENA experiment ERX168813 (http://www.ebi.ac.uk/ena/data/view/ERX168813).
