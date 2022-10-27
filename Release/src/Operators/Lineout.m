function dataStruct = Lineout(dataStruct, coords, meshtype)
%
% Lineout Retrieve a lineout from a two- or three-dimensional dataStruct
%
%-------------------------------------------------------------------------------%
% Info:
%
% Inputs:
%   dataStruct - structure containing all data fields (solution field, mesh data, etc.)
%   meshtype   - desired discretization; uniform or non-uniform
%
% Outputs:
%
% Licensing:
%   This file is part of FLASH-AVM.
%   
%   FLASH-AVM is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   FLASH-AVM is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with FLASH-AVM.  If not, see <https://www.gnu.org/licenses/>.
%
%-------------------------------------------------------------------------------%

  % give error if not 2D or 3D data
  if dataStruct.ndim < 2
    error('Input data dimension must be greater than one.');
  end

  % give error if dimensions of coords incorrect
  if dataStruct.ndim == 2
    if size(coords) ~= [2 2]
      error('Incorrect input dimensions of coords argument. Dimensions should be [2 2].');
    end
  end
  if dataStruct.ndim == 3
    if size(coords) ~= [2 3]
      error('Incorrect input dimensions of coords argument. Dimensions should be [2 2].');
    end
  end

  % minimum cell size in the mesh
  dr = min(dataStruct.nonuniform.meshres) / 2;

  % parameterize lineout
  delx = (coords(2,1) - coords(1,1));
  dely = (coords(2,2) - coords(1,2));
  delr = sqrt(delx^2 + dely^2);
  theta = acos(delx/delr);
  
  % compute discrete lineout array
  lineout_r = 0:dr:delr;
  lineout_x = coords(1,1) + lineout_r * cos(theta);
  lineout_y = coords(1,2) + lineout_r * sin(theta);

  % initialize lineout cell tracker 
  lineout_cells = false(length(dataStruct.nonuniform.mesh(:,1,1))-1, ...
                        length(dataStruct.nonuniform.mesh(:,2,1))-1, ...
                        1, ...
                        dataStruct.nonuniform.nblocks);

  % set default discretization style
  if nargin < 2
    meshtype = 'nonuniform';
  end

  % get number of lineouts currently in the datastructure
  if isfield(dataStruct, 'lineout')
    nlo = length(dataStruct.lineout);
  else
    nlo = 0;
  end

  % produce 1D curve based on desired meshtype
  if strcmp(meshtype, 'nonuniform')

    % initialize curve arrays
    dataStruct.lineout{nlo+1}.r = [];
    dataStruct.lineout{nlo+1}.xc = [];
    dataStruct.lineout{nlo+1}.yc = [];
    dataStruct.lineout{nlo+1}.dx = [];
    for v = 1:dataStruct.nvars
      dataStruct.lineout{nlo+1}.(dataStruct.vars{v}) = [];
    end

    % loop through elements of discrete lineout array
    for m = 1:length(lineout_r)

      % loop through blocks in data structure
      for blk = 1:dataStruct.(meshtype).nblocks

        % first check if element is within current block
        if lineout_x(m) >= dataStruct.nonuniform.meshbnd(1,1,blk) ...
          & lineout_x(m) <= dataStruct.nonuniform.meshbnd(2,1,blk) ...
          & lineout_y(m) >= dataStruct.nonuniform.meshbnd(1,2,blk) ...
          & lineout_y(m) <= dataStruct.nonuniform.meshbnd(2,2,blk)

          % loop through the cells within the block
          for i = 1:length(dataStruct.nonuniform.mesh(:,1,blk))-1
            for j = 1:length(dataStruct.nonuniform.mesh(:,2,blk))-1
              for k = 1:1

                % get spatial information for this cell
                xlo = dataStruct.nonuniform.mesh(i,1,blk);
                xhi = dataStruct.nonuniform.mesh(i+1,1,blk);
                ylo = dataStruct.nonuniform.mesh(j,2,blk);
                yhi = dataStruct.nonuniform.mesh(j+1,2,blk);
                xc = 0.5*(xlo + xhi);
                yc = 0.5*(ylo + yhi);
                rc = sqrt((xc-coords(1,1))^2 + (yc-coords(1,2))^2);

                % check if lineout element is contained within this cell and we don't already have cell in array
                if lineout_x(m) > xlo & lineout_x(m) <= xhi ...
                  & lineout_y(m) > ylo & lineout_y(m) <= yhi ...
                  & lineout_cells(i,j,k,blk) == false

                  % add this cell to the lineout array
                  lineout_cells(i,j,k,blk) = true;

                  % set cell midpoints and cell size
                  dataStruct.lineout{nlo+1}.r = [dataStruct.lineout{nlo+1}.r rc];
                  dataStruct.lineout{nlo+1}.xc = [dataStruct.lineout{nlo+1}.xc xc];
                  dataStruct.lineout{nlo+1}.yc = [dataStruct.lineout{nlo+1}.yc yc];
                  dataStruct.lineout{nlo+1}.dx = [dataStruct.lineout{nlo+1}.dx (xhi-xlo)];

                  % set variable data
                  for v = 1:dataStruct.nvars
                    dataStruct.lineout{nlo+1}.(dataStruct.vars{v}) = [dataStruct.lineout{nlo+1}.(dataStruct.vars{v}) ...
                      dataStruct.nonuniform.(dataStruct.vars{v})(i,j,k,blk)];
                  end

                end

              end
            end
          end

        end % end block check loop

      end % end line element loop

    end

  else

  end

end
