function [] = PlotMesh2D(dataStruct,color)
%
% PlotMesh2D Plot the mesh block outlines
%
%-------------------------------------------------------------------------------%
% Info:
%
% Inputs:
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
    for k = 1:dataStruct.nonuniform.nblocks

      % plot bottom line
      px = linspace(dataStruct.nonuniform.meshbnd(1,1,k),dataStruct.nonuniform.meshbnd(2,1,k),np);
      py = dataStruct.nonuniform.meshbnd(1,2,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

      % plot top line
      py = dataStruct.nonuniform.meshbnd(2,2,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

      % plot left line
      py = linspace(dataStruct.nonuniform.meshbnd(1,2,k),dataStruct.nonuniform.meshbnd(2,2,k),np);
      px = dataStruct.nonuniform.meshbnd(1,1,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

      % plot right line
      px = dataStruct.nonuniform.meshbnd(2,1,k) * ones(1,np);
      plot3(px,py,1e30*ones(1,np),'color',color); hold on;

    end
    hold off;

end
