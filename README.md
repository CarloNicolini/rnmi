# rnmi
Relative Normalized Mutual Information wrappers for Matlab/Octave

This folder include c++ code to compute relatieve Normalized Mutual Information of two partitions. rNMI is Normalized Mutual Information (NMI) of given two partitions minuses the expected NMI of two random partitions, with the same number of groups as given partitions. rNMI overcomes the finite-size effect of NMI. You can find more information in this paepr: http://arxiv.org/abs/1501.03844

The usage is easy: after compiling using “make”, there will be 2 executable files “nmi” and “rnmi”, which calculate NMI and rNMI respectively. 
to calculate NMI of two configuration stored in 1.conf and 2.conf (labeling start form 0), you can use
“./nmi 1.conf 2.conf”

To compute rNMI, 
“./rnmi 1.conf 2.conf”

Pan Zhang
pan@santafe.edu
http://panzhang.net

Carlo Nicolini
