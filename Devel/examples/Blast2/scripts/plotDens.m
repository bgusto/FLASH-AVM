clc, clear all, clf

% add flash-avm
addpath(genpath(getenv('FLASHAVM_DEV')));

% set latex
set(0,'defaulttextinterpreter','latex');

% plotting options
fs = 16;

% load data
dataStruct = GrabData('../data/blast2_hdf5_chk_0018', {'dens'});

% apply 'curve' operator to create 1D arrays
dataStruct = Curve(dataStruct, 'nonuniform');

% plot curve
plot(dataStruct.curve.xc, dataStruct.curve.dens, 'k.');
xlabel('$x$', 'fontsize', fs);
ylabel('$\rho$', 'fontsize', fs);
