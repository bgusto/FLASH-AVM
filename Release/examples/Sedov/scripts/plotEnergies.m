clc, clear all, clf

% add flash-avm to path
addpath(getenv('FLASHAVM'));

% misc. options
fs_labl = 16;       % abcissa-ordinate label fontsize
fs_axes = 10;       % fontsize for tick labels
fs_lgnd = 12;       % fontsize for legend
lw = 1.5;

% name of integrated quantities file
filenm = '../data/sedov.dat';

% get cell count from integral quantities file
iq = GrabIntegralQuantity(filenm);

% plot the number of cells over time
semilogy(iq.time,iq.E_kinetic,'b','linewidth',lw); hold on;
semilogy(iq.time,iq.E_internal,'r','linewidth',lw); hold on;
semilogy(iq.time,iq.E_total,'k','linewidth',lw);

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
