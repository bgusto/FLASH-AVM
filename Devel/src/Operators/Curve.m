function dataStruct = Curve(dataStruct, meshType, refLevel)
%
% Curve Convert a dataStruct from 'GrabData' to a one-dimensional curve
%
%-------------------------------------------------------------------------------%
% Info:
%   This function takes the (one-dimensional) data output from GrabData.m
%   and reformats the nonuniform (multi-block) data into a contiguous array. The
%   output data will be nonuniform if the original data set uses adaptive mesh
%   refinement unless the input field 'refLevel' is set to the desired output AMR
%   level, in which case the nonuniform data will be prolonged to that desired
%   level. The routine creates a new field in the 'dataStruct' structure called
%   'curve'.
%
% Inputs:
%   dataStruct - structure containing all data fields (solution field, mesh data, etc.)
%   meshType   - desired discretization; uniform or non-uniform
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

  % give error if not 1D data
  if dataStruct.ndim ~= 1
    error('Input data must be one-dimensional.');
  end

  % set default discretization style
  if nargin < 2
    meshType = 'nonuniform';
  end

  % produce 1D curve based on desired meshType
  if strcmp(meshType, 'nonuniform')

    % initialize curve arrays
    dataStruct.curve.x = [];
    dataStruct.curve.xc = [];
    dataStruct.curve.dx = [];
    for v = 1:dataStruct.nvars
      dataStruct.curve.(dataStruct.vars{v}) = [];
    end

    % order the blocks from xmin to xmax
    [tmp, ord] = sort(dataStruct.nonuniform.meshbnd(1,1,:));

    % loop through blocks in data structure
    for b = 1:dataStruct.nonuniform.nblocks

      % get block id from ordered array
      blk = ord(b);

      % set cell edges and midpoints
      dataStruct.curve.x = [dataStruct.curve.x dataStruct.nonuniform.mesh(:,1,blk)'];
      dataStruct.curve.xc = [dataStruct.curve.xc 0.5*(dataStruct.nonuniform.mesh(1:end-1,1,blk)' ...
        + dataStruct.nonuniform.mesh(2:end,1,blk)')];

      % set mesh spacing
      dataStruct.curve.dx = [dataStruct.curve.dx ...
                              dataStruct.nonuniform.meshres(blk) ...
                              * ones(1,length(dataStruct.nonuniform.mesh(:,1,blk)))];

      % set variable data
      for v = 1:dataStruct.nvars
        dataStruct.curve.(dataStruct.vars{v}) = [dataStruct.curve.(dataStruct.vars{v}) ...
          dataStruct.nonuniform.(dataStruct.vars{v})(:,1,1,blk)];
      end

    end

    % finally sort data based on coordinates
    dataStruct.curve.x = unique(dataStruct.curve.x);
    [dataStruct.curve.xc ord] = sort(dataStruct.curve.xc);
    dataStruct.curve.dx = dataStruct.curve.dx(ord);
    for v = 1:dataStruct.nvars
      dataStruct.curve.(dataStruct.vars{v}) = dataStruct.curve.(dataStruct.vars{v})(ord);
    end

  elseif strcmp(meshType, 'uniform')

    % get max mesh level
    if nargin < 3
      lref_max = max(dataStruct.nonuniform.lrefine);
    else
      lref_max = refLevel;
    end

    % first we have to count the size of required arrays :(
    nxb = length(dataStruct.nonuniform.mesh(:,1,1)) - 1;
    cnt = 0;
    for b = 1:dataStruct.nonuniform.nblocks
      lref = dataStruct.nonuniform.lrefine(b);
      ldiff = lref_max - lref;
      cnt = cnt + 2^ldiff * nxb;
    end 

    % initialize curve arrays
    dataStruct.curve.x = zeros(1,cnt);
    dataStruct.curve.xc = zeros(1,cnt);
    dataStruct.curve.dx = zeros(1,cnt);
    for v = 1:dataStruct.nvars
      dataStruct.curve.(dataStruct.vars{v}) = zeros(1,cnt);
    end

    % order the blocks from xmin to xmax
    [tmp, ord] = sort(dataStruct.nonuniform.meshbnd(1,1,:));

    % loop through blocks in data structure
    cnt = 0;
    for b = 1:dataStruct.nonuniform.nblocks

      % get block id from ordered array
      blk = ord(b);

      % get the mesh level of the block
      lref = dataStruct.nonuniform.lrefine(blk);
      ldiff = lref_max - lref;

      % geom factors
      xlo = dataStruct.nonuniform.meshbnd(1,1,blk);
      xhi = dataStruct.nonuniform.meshbnd(2,1,blk);
      ncells = (2^ldiff) * length(dataStruct.nonuniform.mesh(1:end-1,1,blk));
      xloc = linspace(xlo, xhi, ncells+1);

      % set cell edges and midpoints
      dataStruct.curve.x(cnt+(1:ncells+1)) = xloc;
      dataStruct.curve.xc(cnt+(1:ncells)) = 0.5*(xloc(1:end-1) + xloc(2:end));

      % set mesh spacing
      dataStruct.curve.dx(cnt+(1:ncells)) = dataStruct.nonuniform.meshres(blk) * ones(1,ncells);

      % loop through vars
      for v = 1:dataStruct.nvars

        % loop through original cells
        for i = 1:nxb;

          ilo = 2^ldiff * (i-1) + 1;
          ihi = 2^ldiff * i;
          dataStruct.curve.(dataStruct.vars{v})(cnt+(ilo:ihi)) = ones(1,2^ldiff) * dataStruct.nonuniform.(dataStruct.vars{v})(i,1,1,blk);

        end

      end

      cnt = cnt + ncells;

    end

    %% finally sort data based on coordinates
    %dataStruct.curve.x = unique(dataStruct.curve.x);
    %[dataStruct.curve.xc ord] = sort(dataStruct.curve.xc);
    %dataStruct.curve.dx = dataStruct.curve.dx(ord);
    %for v = 1:dataStruct.nvars
    %  dataStruct.curve.(dataStruct.vars{v}) = dataStruct.curve.(dataStruct.vars{v})(ord);
    %end

  end

end
