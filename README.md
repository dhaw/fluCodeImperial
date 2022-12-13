# fluCodeImperial
All codes used for "Using real-time data to guide decision-making during an influenza pandemic: a modelling analysis" are included int his folder. The file runExamples.m contains a step-by-step method for reproducing figures and running model fits. Ensure that all data and code files are in the same directory, then a single execution of "runExamples" on the command line will generate all main and supplementary figures in the manuscript. There is one line of code per figure, clearly marked, that can be commented out as desired. In order to run the MCMC adaptive algorith, a single line of code, also clearly marked, must be commented back in. Instructions to change the single-state example are given at the top of the file "runExamples.m". Selecting California ("state=1") is consistent with allr esults presented in the manuscript. 

Plots make use of files from the following sources, with some modifications:
Holger Hoffmann (2022). Violin Plot (https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot);
Evan (2022). Plot Groups of Stacked Bars (https://www.mathworks.com/matlabcentral/fileexchange/32884-plot-groups-of-stacked-bars);
John Onofrey (2022). Shaded Plots and Statistical Distribution Visualizations (https://www.mathworks.com/matlabcentral/fileexchange/69203-shaded-plots-and-statistical-distribution-visualizations)
