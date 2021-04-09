function [] = PlotMesh2D(mesh,color)
%
% PlotMesh2D Plot the mesh block outlines
%
%-------------------------------------------------------------------------------%
% Info:
%
% Inputs:
%   mesh  - the array containing mesh boundary info
%   color - desired color, specified as rgb string (i.e. '[0.1 0.1 0.1]') for
%           block outlines (default is gray)
%
% Outputs:
%
%-------------------------------------------------------------------------------%

    % optional arguments (default color)
    if nargin < 2
        color = '[0.8 0.8 0.8]';
    end 

    % number of points per line
    np = 2;

    % plot mesh
    for k = 1:length(mesh)

      % width, length, and number of points to use
      w = mesh(2,k) - mesh(1,k);
      h = mesh(4,k) - mesh(3,k);

      % plot bottom line
      px = linspace(mesh(1,k),mesh(2,k),np);
      py = mesh(3,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

      % plot top line
      py = mesh(4,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

      % plot left line
      py = linspace(mesh(3,k),mesh(4,k),np);
      px = mesh(1,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

      % plot right line
      px = mesh(2,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

    end
    hold off;

end
