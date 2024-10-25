clear;clc;close all
addpath(['..', filesep, 'DataDerived'])
T = readtable('20210920_dev9.csv');

% rename variables to match cells format
T = renamevars(T,["x","y","C1","C2","C3"], ...
                 ["X","Y","c1","c2","c3"]);
             
cells = ternaryZoning(T);
rmpath(['..', filesep, 'DataDerived'])