clear;clc;close all

data = readtable('none-MRX-VCR.csv');
% c1 = none
% c2 = MRX
% c3 = VCR
ternaryZoning(data)
figure(2)
legend({'MRX-VCR', 'Media-VCR', 'Media-MRX'})


figure(3)
terlabel('Media', 'MRX', 'VCR')