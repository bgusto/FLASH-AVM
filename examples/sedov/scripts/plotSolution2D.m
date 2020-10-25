clc, clear all, clf

% add flash-avm to path
addpath(getenv('FLASHAVM'));

% set latex
set(0,'defaulttextinterpreter','latex');

% load the hot_desat colormap
imap = load('hot_desaturated.mat');
imap = imap.hot_desat;

% misc. options
fs_labl = 15;       % abcissa-ordinate label fontsize
mt = 0.5;           % mesh line thickness

% set filename and basename
chk = 383;

% file basename
basenm = '../data/sedov_hdf5_plt_cnt_';

% get data
[pres m x y] = GrabData2D([basenm sprintf('%0.4d',chk)],'pres','avginterp');

% viewing window
pxmin = 0.5;
pxmax = 1;
pymin = 0.5;
pymax = 1;

% figure handle
figure(1)

% get axes ranges
xmin = min(x(:));
xmax = max(x(:));
ymin = min(y(:));
ymax = max(y(:));

% number of cells
[ny nx] = size(pres);

% new meshgrid
[xp yp] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

% plot pseudocolor
surf(xp,yp,pres,'edgecolor','none'); view(2); hold on;

% axes
axis([pxmin pxmax pymin pymax]);
%xticks([0.5 0.55 0.6 0.65 0.7]);
%yticks([0.5 0.55 0.6 0.65 0.7]);

% labels
xlabel('$x$ [cm]','fontsize',fs_labl,'interpreter','latex');
ylabel('$y$ [cm]','fontsize',fs_labl,'interpreter','latex');

% options
box on;
colormap(imap)
c = colorbar;
%c.YTick = [3e09:0.75e09:6e09];
%caxis([3.e09 6.e09]);
%caxis([3.e09 6.e09]);
pbaspect([(pxmax-pxmin)/(pymax-pymin) 1 1]);
set(gca, 'Layer', 'top');
grid off

% plot mesh
for k = 1:length(m)

  % width, length, and number of points to use
  w = m(2,k) - m(1,k);
  h = m(4,k) - m(3,k);
  np = 128;

  % plot bottom line
  px = linspace(m(1,k),m(2,k),np);
  py = m(3,k) * ones(1,np);
  plot3(px,py,1e10*ones(1,np),'color','[0.8,0.8,0.8]'); hold on;

  % plot top line
  py = m(4,k) * ones(1,np);
  plot3(px,py,1e10*ones(1,np),'color','[0.8,0.8,0.8]'); hold on;

  % plot left line
  py = linspace(m(3,k),m(4,k),np);
  px = m(1,k) * ones(1,np);
  plot3(px,py,1e10*ones(1,np),'color','[0.8,0.8,0.8]'); hold on;

  % plot right line
  px = m(2,k) * ones(1,np);
  plot3(px,py,1e10*ones(1,np),'color','[0.8,0.8,0.8]'); hold on;

end
hold off;

% print figure
%print(figure(1),'sedov-n4096-amr-eps1e01-pressure-mesh','-dpng','-r700');
