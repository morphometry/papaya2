% You may need to modify the "mex" line to add include and library paths.
% see imt_for_pointpattern.cpp for more information
mex -v CXXFLAGS='-std=c++14' -I"../include" -I"/usr/local/include" -L"/usr/local/lib" -lgmp imt_for_pointpattern.cpp;
a = readtable('../demos/example_inputs/granular-cryst-cluster.txt');
a = table2array(a);

% without periodic boundary conditions
imt_open = imt_for_pointpattern(a)
% with periodic boundary conditions
imt_pbc = imt_for_pointpattern(a, [500., 500.])
% look at the q6 values for all the seeds
imt_pbc(:,9)

figure;
subplot(1,2,1);
voronoi(a(:,1),a(:,2));
subplot(1,2,2);
hist(imt_pbc(:,9),30);
xlabel('q_6');
ylabel('abundance');
