DiffLogo
========
What is DiffLogo
---------------
The DiffLogo tool is a R package for the visualization of differences between multiple motifs for different alphabets. The user supplies a set of motifs each represented as position weight matrices (PWMs) [1]. The DiffLogo package supports the comparison of two motifs by a single DiffLogo and the comparison of multiple motifs by a DiffLogo-table. DiffLogo is based on the idea behind the well-known sequence logo [2], i.e. motifs are visualized position-wise based on two functions. First, the <code>stackHeight</code> function computes the height of each stack. Second, the <code>baseDistribution</code> function breaks down the stack height on the individual characters. The user is able to parametrise the individual functions with arbitrary functions <code>stackHeight</code> and <code>baseDistribution</code>. Default implementations of both functions are provided.

Who should use DiffLogo
-----------------------
DiffLogo is designed for researchers with any biological background and interest in computational biology. It is useful to document findings, share knowledge, and to present the outcome of motif prediction pipelines. It aims at an increase of the quality of comparative publications. In addition, DiffLogo eases the communication between bioinformatics scientists and biologists.

When should DiffLogo be used
----------------------------
DiffLogo is intended for the comparison of similar motifs. These motifs can come from different sources such as
- different treatments, species, cell lines, or tissues
- different motif prediction algorithms and configurations

Getting started
---------------
Download and install the R package DiffLogo available in the folder 'built'. Load your motifs of interest as PWMs [1] to R. Please find example motifs in the file 'exampleLogos.R' and in the folder 'inst/pwm' (extracted from [3]). Please find the vignette 'DiffLogoBasics.Rnw' of the DiffLogo package for example code.  

Use the function <code>diffLogoFromPwm(PWM1, PWM2)</code> to visualize the difference between two motifs. 

Use the function <code>diffLogoTable(listOfPWMs)</code> to visualize the pairwise difference between more than two motifs.

References
--------------

[1] http://en.wikipedia.org/wiki/Position_weight_matrix<br>
[2] Schneider TD, Stephens RM. 1990. Sequence Logos: A New Way to Display Consensus Sequences. _Nucleic Acids Res. 18_:6097-6100<br>
[3] Eggeling, R., Gohr, A., Keilwagen, J., Mohr, M., Posch, S., Smith, A.D., Grosse, I.: On the value of intra-motifdependencies of human insulator protein ctcf. _PLoS ONE 9(1)_, 85629 (2014). doi:10.1371/journal.pone.0085629
