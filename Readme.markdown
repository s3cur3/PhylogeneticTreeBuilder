Instructions for using the Phylogenetic Tree Builder
==========================================

To run the phylogenetic tree builder, first make sure you're in a *nix
environment (Linux or Mac OS). Then, do the following:

1. Unzip the zip containing this program's source.
2. Obtain some FASTA file whose sequences you wish to use to build a phylogenetic tree.
3. Run `team_2_main.py`, supplying it the name of your FASTA file and the number of bootstrapping iterations you wish to run.
    - `$ python3 team_2_main.py --bootstrap <num bootstrapping iterations> <your FASTA file name>.fasta`
	- If you have problems, get help by running:
	    - `$ python3 team_2_main.py -h`
4. The program will churn for a long time, first to generate all pairwise sequence comparisons (which it will write to a file in your working directory, named `<your FASTA file name>.txt`) and then to bootstrap your aligments.
5. When it finishes, the program will tally up the number of times that a given phylogenetic tree occurred. To visualize the tree, copy its text (which looks something like `(H. sapiens Navaho, (H. sapiens German, H. sapiens Japan));`) into a `.tre` file and open it in a tree visualization program like [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). 

