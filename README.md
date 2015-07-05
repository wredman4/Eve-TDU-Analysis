# Eve-TDU-Analysis
# Code for analysis done on Eve looking at topological distance unit (TDU) relationships. 

This repository is a colloboration between William Redman (NYU & Princeton), Danny Wells (Northwestern & Princeton), and Madhav Mani.
Please note that this repository includes the SLM toolbox developed by John D'Errico (http://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling)
in its unabrigded form (including license), which is the only condition for redistribution. Please see SLM documentation for more information.

Eve TDU Analysis (ETA) has five components:

1) Mutual Information Analysis (MIA): MutualInformation_v3.m and MI_v3.m

2) Mutual Information Movie Maker (MIMM): EveMIMovieMaker.m and MutualInfoMovie.m

3) Nuclei Turn-Off Analysis (NTA): TurnOff_v3.m

4) Preprocessing Data Structure Initialization (PDSI): ReadGregorEvePipeline.m

5) Wave Analysis (WA): EveWaveAnalysis.m

The PDSI takes in the data structures created by the FISH Toolbox and LivemRNA and creates some new data structures, namely relationships between different nuclei and particles.
The MIA was done after the integrated, detrended correlation analysis that worked well on the Hunchback data failed to go anywhere.
The MIMM gives both a very visual representation of what is going on through time, as well as helping determine the time scale of the interactions.
The NIA looks at how closely neighboring nuclei's turn-off "profile" is.
The WA looks at both the width of the stripe through time, and the velocity of the stripe. 

Any questions regarding this repository can be sent to wtr214@nyu.edu. 


