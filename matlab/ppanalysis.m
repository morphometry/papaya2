% add -I (includepath) -L (libpath) -l (libraries) as required
% see imt_for_pointpattern.cpp for more information
mex -v -lCGAL -lCGAL_Core -lgmp imt_for_pointpattern.cpp
a = readtable('test_pp.txt');
a = table2array(a);
% without periodic boundary conditions
imt_open = imt_for_pointpattern(a)
% with periodic boundary conditions
imt_pbc = imt_for_pointpattern(a, [500., 500.])
% look at the q2 values for all the seeds
imt_pbc(:,5)
