clc, clear all, clf

% add flash-avm to path
addpath(getenv('FLASHAVM'));

% misc. options
fs_labl = 18;       % abcissa-ordinate label fontsize
fs_axes = 11;       % fontsize for tick labels
fs_lgnd = 14;       % fontsize for legend
lw = 1.5;

% name of integrated quantities file
filenm = '../data/sedov.dat';

% get cell count from integral quantities file
[ekin time] = GrabIntegralQuantity(filenm,'E_kinetic');
[eint time] = GrabIntegralQuantity(filenm,'E_internal');
[etot time] = GrabIntegralQuantity(filenm,'E_total');

% plot the number of cells over time
semilogy(time,ekin,'b','linewidth',lw); hold on;
semilogy(time,eint,'r','linewidth',lw); hold on;
semilogy(time,etot,'k','linewidth',lw);

% options

% axes
axis([-0.002 0.102 1e-1 2]);
set(gca,'fontsize',fs_axes);

% labels
xlabel('$t$','fontsize',fs_labl,'interpreter','latex');
ylabel('Energy','fontsize',fs_labl,'interpreter','latex');

% legend
legend({'Kinetic','Internal','Total'}, ...
    'location', 'northeast', ...
    'interpreter','latex', ...
    'fontsize', fs_lgnd);
legend boxoff;
