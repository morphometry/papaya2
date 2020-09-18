% You may need to modify the "mex" line to add include and library paths.
% see imt_for_pointpattern.cpp for more information
mex -v CXXFLAGS='-std=c++14' -I"../include" -lgmp imt_for_pointpattern.cpp;
a = readtable('test_pp.txt');
a = table2array(a);
% without periodic boundary conditions
imt_open = imt_for_pointpattern(a)
% with periodic boundary conditions
imt_pbc = imt_for_pointpattern(a, [500., 500.])
% look at the q2 values for all the seeds
imt_pbc(:,5)
