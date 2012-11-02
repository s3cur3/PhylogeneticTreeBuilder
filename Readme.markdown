Instructions for Running the Bootstrap Alignment Tool
====================================================

To create a set of bootstrap (randomized) sequences, get all pairwise alignments, and write the data to your hard drive, do the following:

1. Download and unzip the ZIP file containing this project's code. (Use the link near the top of this page.)
2. Open a command prompt (a.k.a. a Terminal window)
3. Type the following (where the `$` indicates the command line waiting for input):
	- `$ cd <whatever directory you unzipped the files to>`
	- `$ python3 bootstrap_and_align.py mtDNA.fasta`
4. The program will read in all the human mitochondrial DNA sequences and create new random sequences therefrom. Then, it will churn for a long time (on my laptop, I estimate it will take about 7 minutes per comparison, or 14 hours for all 120 comparisons). When it finishes, you should have two files, called `bootstrap<some number>.fasta` and `alignments<some number>.txt` in your directory. Email those to Tyler, and go back to step 2!

By my estimates, the program will require a little over 1 GB of RAM when running. Depending on the amount of RAM in your computer, you can potentially run many instances of the program at once (thus giving us more data in a shorter amount of time). The following will tell you how many instances you can run at once, depending on the amount of RAM in your computer:

- 2 GB: 1 instance at a time
- 3 GB: 2 instances at a time
- 4 GB: 3 instances at a time

To run many instances of the program at once, just open up a new command prompt/terminal window and continue from Step 3 above.

If you absolutely cannot tolerate your computer being a bit sluggish, just run that number of instances overnight (or subtract 1 from that number). I advise that you do *not* allow your computer to sleep or hibernate while running the computer---I'm not sure how it will behave.

To do neighbor joining properly, we'd like as many bootstrapping runs as possible. If you can contribute just 1 run per day, we can get 50 data points when writing our paper, which would be *awesome*!


Instructions for using the Phylogenetic Tree Builder
==========================================
**At this point, this is way not ready to use. Check back later.**

To run the drift simulator, first make sure you're in a *nix
environment (Linux or Mac OS). Then, do the following:

1. Unzip the zip.
2. . . .
3. Profit?

