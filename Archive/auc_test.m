% script to test the gradient AUC method against other synergy
% based on a full dose matrix
clear;clc

data = readtable('testdrugs.csv');
[c1_mat, c2_mat, doseMatrix] = drugResponseVec2Mat(data.C1, data.C2, data.inhibition);

% rows going down are doses in concentration 1
% columns going to right are doses in concentration 2

diagonal_dose = diag(doseMatrix);