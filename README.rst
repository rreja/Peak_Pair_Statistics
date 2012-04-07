<<<<<<< HEAD
Peak Statistics
=======
Peak-pair Statistics
>>>>>>> d07c84c6eb50aec923beb411a93c3b718abdfbe3
================

Introduction
-------------

<<<<<<< HEAD
The scriptis in this repository can be used to perform some basic statistics
on the "peak calls". The "peak calls" file is the output file obtained by running
genetrack program on the index/raw tag file.

The basic operations include seperating the singleton and non-singelton peaks and
then calculating the median and average of tag counts and standard deviations.

=======
The scripts in this repository can be used to perform some basic statistics on the peak-pairs. The peak-pair file is the output file obtained by running cw-peak-pair python script on the peak call file.

The basic operations include calculating a) the mean and median of peak-pair occupancy, b) peak-pair mode, c) Counting the number of orphans, d) Fraction of all mapped reads that reside in peak-pairs and, e) Signal to noise ratio in the dataset.
>>>>>>> d07c84c6eb50aec923beb411a93c3b718abdfbe3

Requirements
------------

<<<<<<< HEAD
- The software requires only Perl_ (5 or higher) to run.
- The input file should be in standard Gff_ format

Installation
------------

Unpack the source code archive. The folder contains the followng

- `peak_statistics.pl`: Script for basic statistics when running on a single file
- `peak_statistics_batch.sh` : Script for batch processing
- `README.rst` : Readme file
- `sample.gff`: The sample input file to test the scripts.

When running the script on a single input gff file:

- Open the file `peak_statistics.pl` in  any text editor of your choice.
- Comment out the line 18 by putting a "#" in front of it.
- Uncomment the lines 21 and 22 by removing the "#" in front of them.
- You are ready to use the file. How to run the script from your terminal?

- Type the following:

    $ perl  peak_statistics.pl  <path_to_your_input_file>
    $ #for example: peak_statistics.pl /Users/input/sample.gff

An output file with  '_NoS.gff' at the end will be generated, that will
contain all the non-singelton peaks into it (The peaks whose standard deviation
is greater than 0). Also, some statistic summary about the file would be displayed
on the screen. Each column (in the order of display)  contains the following:
    - Filename
    - Percentage of mapped reads
    - Percentage of uniquely mapped reads
    - Total non-singelton peaks
    - Total singelton peaks
    - Median of tag counts for non-singleton peaks
    - Mean of tag counts for non-singleton peaks
    - Median fuzziness (standard deviation) for non-singleton peaks
    - Mean fuzziness (standard deviation) for non-singleton peaks

When batch processing:

- For batch processing DO NOT modify 'peak_statistics.pl' file
- To run the batch script:

    $./peak_statistics_batch.sh # will give the instruction on how to run the file
    $./peak_statistics_batch.sh -i <path_to_directory>
    $ # for example: ./peak_statistics_batch.sh  /Users/input/

Remember, when giving the path to the directory containing files, always give the
path with "/" at the end. For ex: "/Users/input/" is acceptable but "/Users/input"
would generate error.


The batch script would create a foler "output" in the directory containing your
input files. The output folder will contain a "_NoS.gff" file for each input
gff file and a "peak_statistics.txt" file containing the basic statistics for
all the input files present in the directory. "peak_statistics.txt" has the
same format as mentioned above.


.. _Perl: http://www.perl.org/
.. _Gff: http://genome.ucsc.edu/FAQ/FAQformat#format3
~                                                                                                                                                                                                                                                                             

=======
- The script only requires Perl_ (5 or higher) to run.
- The input tag file should have the idx/tab extension and should be of the form (chr,index,forward,reverse,value(optional column)).
- The peak-pair files should be in the standard gff format.


THE SCRIPT WILL BREAK IF:
------------------------

- The files have excel ^M character in it. For sanity check, open your file in terminal, to see if you can see ^M character in your file. In case, you find ^M character in your file, use the following command to remove it::

    $ perl -p -e 's/^M/\n/g;' <file_with_excel_char> > <new_file>

- The peak-pair file should start with a "S_" and should end with a "_sXeXFX.gff", where X could be any number, ex: "_s5e20F10.gff".
- The orphan file should start with a "O_".
- The names of all S_* and O_* files should contain the index file name in it. For ex. if index file name is "Reb1-rep2.idx", than all the S_* and O_* files should be like S_XXX_Reb1-rep2_XXX_sXe20F1.gff, where X is any character.


Installing and Running the scripts
------------

Unpack the source code archive. The folder contains the following::

-  robust_peak_pair_stats.pl: Script for basic statistics and an increasing-window quantile scan for signal:noise. 
-  pp_stats_5pt_scan.pl: Script for basic statistics and a fixed-window quantile scan for signal:noise.
-  README.rst: Readme file
-  Sample data: which includes (two index files: Reb1-rep2.idx and Reb1-rep3.idx) and folder (genetrack_s5e10F1) containing peak calls and a subfolder (cwpair_output_mode_f0u0d100b3) containing all the S_*, D_*, O_*, and P_*, peak-pair files


If you wish to get the signal:noise ratio infomation using increasing-window quantile scan (for ex. top 1%, top 5%, top 10% etc) than use the following script::

    $ perl  robust_peak_pair_stats.pl -h

Do a test run of the script by typing::

$ perl robust_peak_pair_stats.pl -i  ./ -d genetrack_s5e10F1/cwpair_output_mode_f0u0d100b3/ -g sg07

The folder should now contain a "peak_pair_stats.txt" file. This means that script runs fine on your system.

if you wish to get the signal:noise information using fixed-width quantile scan (for ex. 0-5 %, 5-10 %, 10- 15 %) than use the following script::

    $ perl pp_stats_5pt_scan.pl -h

Do a test run of the script by typing::

    $  perl pp_stats_5pt_scan.pl -i  ./ -d genetrack_s5e10F1/cwpair_output_mode_f0u0d100b3/ -g NA -s 160000000 -p 10

The folder should now contain, a "peak_pair_stats.txt" and a "signal2noise_qt_scan.txt" file.
This means that script runs fine on your system.


Output
------

All output files will be produced in the folder that contain S_* and O_* files.
Following output files will be generated:

- The script "pp_stats_5pt_scan.pl" produces an extra  file named: "signal2noise_qt_scan.txt", which will contain the quantile range and the signal to noise ratio in a tab delimited format.

- "peak_pair_stats.txt" containing the summary for each input file. The summary includes the following information::

    - Filename
    - Peak-pair mode
    - Peaks in peak pairs
    - Orphan peaks
    - Median peak-pair occupancy
    - Mean peak-pair occupancy
    - FRIP (Fraction of all mapped reads in peak-pairs) 
    - top_1pt_signal:noise [only in the output of "robust_peak_pair_stats.pl"]
    - top_5pt_signal:noise [only in the output of "robust_peak_pair_stats.pl"]
    - top_10pt_signal:noise [only in the output of "robust_peak_pair_stats.pl"]
    - top_25pt_signal:noise [only in the output of "robust_peak_pair_stats.pl"]
    - top_50pt_signal:noise [only in the output of "robust_peak_pair_stats.pl"]
    - top_75pt_signal:noise [only in the output of "robust_peak_pair_stats.pl"]
    - top_100pt_signal:noise [only in the output of "robust_peak_pair_stats.pl"]
	
   

.. _Perl: http://www.perl.org/
.. _Gff: http://genome.ucsc.edu/FAQ/FAQformat#format3
>>>>>>> d07c84c6eb50aec923beb411a93c3b718abdfbe3
