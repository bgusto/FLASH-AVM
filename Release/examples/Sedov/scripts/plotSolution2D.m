clc, clear all, clf

% add flash-avm to path
addpath(genpath(getenv('FLASHAVM')));

% set latex
set(0,'defaulttextinterpreter','latex');

% load the hot_desat colormap
imap = load('hot_desaturated.mat');
imap = imap.hot_desat;

% misc. options
fs_labl = 17;       % abcissa-ordinate label fontsize
mt = 0.5;           % mesh line thickness

% set filename and basename
chk = 383;

% file basename
basenm = '../data/sedov_hdf5_plt_cnt_';

% get data
solnData = GrabData([basenm sprintf('%0.4d',chk)],{'dens','pres'},'copy');

%%%% viewing window
%%%pxmin = 0.5;
%%%pxmax = 1;
%%%pymin = 0.5;
%%%pymax = 1;
%%%
%%%% figure handle
%%%figure(1)
%%%
%%%% get axes ranges
%%%xmin = min(x(:));
%%%xmax = max(x(:));
%%%ymin = min(y(:));
%%%ymax = max(y(:));
%%%
%%%% number of cells
%%%[ny nx] = size(pres);
%%%
%%%% new meshgrid
%%%[xp yp] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
%%%
%%%% plot pseudocolor
%%%surf(xp,yp,pres,'edgecolor','none'); view(2); hold on;
%%%
%%%% axes
%%%axis([pxmin pxmax pymin pymax]);
%%%%xticks([0.5 0.55 0.6 0.65 0.7]);
%%%%yticks([0.5 0.55 0.6 0.65 0.7]);
%%%
%%%% labels
%%%xlabel('$x$ [cm]','fontsize',fs_labl,'interpreter','latex');
%%%ylabel('$y$ [cm]','fontsize',fs_labl,'interpreter','latex');
%%%
%%%% options
%%%box on;
%%%colormap(imap)
%%%c = colorbar;
%%%%c.YTick = [3e09:0.75e09:6e09];
%%%%caxis([3.e09 6.e09]);
%%%%caxis([3.e09 6.e09]);
%%%pbaspect([(pxmax-pxmin)/(pymax-pymin) 1 1]);
%%%set(gca, 'Layer', 'top');
%%%grid off
%%%
%%%% plot mesh
%%%PlotMesh2D(m);
%%%
%%%% print figure
%%%print(figure(1),'sedov-pressure-mesh','-dpng','-r450');
