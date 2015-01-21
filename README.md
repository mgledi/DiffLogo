What is DiffLogo
================
DiffLogo is a tool to visualize differences between two motifs. It is based on the commonly known sequence logo [1]. The only precondition that must been fulfilled in the current version is that the two input motifs must be of the same length. The current version of DiffLogo provides two different visualizations.

Who should use DiffLogo
=======================
DiffLogo is for any researcher with any biological background. It is useful to document your findings, to share your knowledge, and to present the outcome of motif prediction pipelines. It should increase the quality of comarative publications. Not least DiffLogo should ease the communication between bioinformatics scientists and biologists.

When should DiffLogo be used
============================
DiffLogo should be used, when it is needed to compare very similar sequence logos. Those sequence logos can come from e.g. 
- different species, cell lines or, tissues
- different motif prediction algorithms that have been run on the same dataset 
- different wet lab experiments and different preconditions of experiments (e.g. stress)

Getting started
===============
Download the file <a href="diffSeqLogo.R">diffSeqLogo.R</a> to your working directory. Start R and load the file with the command <code>source</code>. Load your motifs of interest as PWMs[2] to R. Please see <a href="exampleLogos.R">exampleLogos.R</a> for some example motifs.  

Use the function <code>diffSeqLogo(PWM1, PWM2)</code> to visualize the difference between two sequence logos. 

Use the function <code>diffSeqLogoMulti(listOfPWMs)</code> to visualize the pairwise difference between more than two motifs.

Examples of usage can be found in the file <a href="compareCTCF.R">compareCTCF.R</a>. The used motifs are provided by Eggeling et al. [3]

Cite this work
==============


[1] Schneider TD, Stephens RM. 1990. Sequence Logos: A New Way to Display Consensus Sequences. _Nucleic Acids Res. 18_:6097-6100

[2] http://en.wikipedia.org/wiki/Position_weight_matrix

[3]
