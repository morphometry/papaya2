% add -I (includepath) -L (libpath) -l (libraries) as required
mex -v -lCGAL -lCGAL_Core -lgmp imt_for_pointpattern.cpp
a = readtable('test_pp.txt');
a = table2array(a);
% without periodic boundary conditions
imt = imt_for_pointpattern(a)
% with periodic boundary conditions
imt = imt_for_pointpattern(a, [500., 500.])
% look at the q2 values for all the seeds
imt(:,5)
