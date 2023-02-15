RMST_Psedo_Evaluation


This project contains the code to produce the results in the manuscript titled, "Estimating Differences in Restricted Mean Survival Time in R with Two New Implementations"

The programs to produce the results are in the 'Programs' folder.

 - the paper_figs.R program produces the output for the figures in the manuscript as well as the example R output. First, run the entire program. Then:

    - figure 1 is produced by lines 90-99. It is saved in the documentation folder with filename RMSTdiffpic.png.
    - the example dataset output at the bottom of page 9 is produced by lines is output by lines 39-42. It is saved in the documentation folder with filename output1.txt.
    - figure 2 is output by lines 56-61.It is saved in the documentation folder with filename km1.png.
    - the R output on page 10 is output by lines 44-53. It is saved in the documentation folder with filename output2.txt.
    - the output on page 11 is output by lines 77-87. It is saved in the documentation folder with filename output4.txt.
    - Figure 3 is produced by lines 104-107. It is saved in the documentation folder with filename kmcurve1.png.
    - Figure 4 s produced by lines 152-188. It is saved in the documentation folder with filename kmcurve2.png.

Tables 1 and 2 are produced by the analysis.R program.


 - run analysis.R. Table 1 is produced by lines 132-205. Table 2 is produced by lines 232-305.


 In order to produce the result datasets, run programs sim1.R, sim1type1.R, sim2.R, sim2type1.R.

 The program sourceRMSTdiff.R is the source file containing the R functions. debug.R is an acillary program for debugging purposes and is not needed for the results. RMST_present.R is a program I used for a presentation and is not necessary for results in the manuscript.
