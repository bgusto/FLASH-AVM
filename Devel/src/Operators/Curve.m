function dataStruct = Curve(dataStruct, meshtype)
%
% Curve Convert a dataStruct from 'Grabata' to a one-dimensional curve
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

  % give error if not 1D data
  if dataStruct.ndim ~= 1
    error('Input data must be one-dimensional.');
  end

  % set default discretization style
  if nargin < 2
    meshtype = 'nonuniform';
  end

  % produce 1D curve based on desired meshtype
  if strcmp(meshtype, 'nonuniform')

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
    for b = 1:dataStruct.(meshtype).nblocks

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

  else

  end

end
