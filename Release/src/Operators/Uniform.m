function dataStruct = Uniform(dataStruct, useProlong, refLevel, interpType, interpOrder, saveNonUni, precision)
%
% Uniform Convert the nonuniform dataStruct to uniformm with possible prolongation
%
%-------------------------------------------------------------------------------%
% Info:
%
% Inputs:
%   dataStruct - structure containing all data fields (solution field, mesh data, etc.)
%
% Outputs:
%   dataStruct.uniform - contains new uniform dataset
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

  % default arguments
  if nargin < 2
    useProlong = false;
  end
  if nargin < 6
    saveNonUni = false;
  end
  if nargin < 7
    precision = 'single';
  end

  % only consider 2d and 3d
  secdim = true;
  terdim = false;
  if dataStruct.ndim == 1
    error('This routine only considers 2D or 3D data. For converting 1D data from nonuniform to uniform, or prolonging 1D data to a finer level, use Curve.m.');
  elseif dataStruct.ndim == 2
  elseif dataStruct.ndim == 3
    terdim = true;
  end

  % first, check the lrefine array
  lref_min = min(dataStruct.nonuniform.lrefine);
  lref_max = max(dataStruct.nonuniform.lrefine);

  % switch based on prolonging or not
  if useProlong

  else

    % get domain boundaries
    xlo = dataStruct.dombnd(1,1);
    xhi = dataStruct.dombnd(2,1);
    ylo = dataStruct.dombnd(1,2);
    yhi = dataStruct.dombnd(2,2);
    zlo = dataStruct.dombnd(1,3);
    zhi = dataStruct.dombnd(2,3);

    % get the number of cells per block
    nxb = length(dataStruct.nonuniform.(dataStruct.vars{1})(:,1,1,1));
    nyb = length(dataStruct.nonuniform.(dataStruct.vars{1})(1,:,1,1));
    nzb = length(dataStruct.nonuniform.(dataStruct.vars{1})(1,1,:,1));

    % get number of blocks per dimension
    nblockx = 0;
    nblocky = 0;
    nblockz = 0;
    for blk = 1:dataStruct.nonuniform.nblocks

      if terdim

        % along x-axis
        if dataStruct.nonuniform.meshbnd(1,2,blk) == ylo ...
          & dataStruct.nonuniform.meshbnd(1,3,blk) == zlo
          nblockx = nblockx + 1;
        end

        % along y-axis
        if dataStruct.nonuniform.meshbnd(1,1,blk) == xlo ...
          & dataStruct.nonuniform.meshbnd(1,3,blk) == zlo
          nblocky = nblocky + 1;
        end

        % along z-axis
        if dataStruct.nonuniform.meshbnd(1,1,blk) == xlo ...
          & dataStruct.nonuniform.meshbnd(1,2,blk) == ylo
          nblockz = nblockz + 1;
        end

      else

        % along x-axis
        if dataStruct.nonuniform.meshbnd(1,2,blk) == ylo
          nblockx = nblockx + 1;
        end

        % along y-axis
        if dataStruct.nonuniform.meshbnd(1,1,blk) == xlo
          nblocky = nblocky + 1;
        end

      end

    end

    % compute the new uniform mesh (assuming constant block sizes for now)
    xe = linspace(dataStruct.dombnd(1,1), ...
                  dataStruct.dombnd(2,1), ...
                  nblockx*nxb+1);
    ye = linspace(dataStruct.dombnd(1,2), ...
                  dataStruct.dombnd(2,2), ...
                  nblocky*nyb+1);
    ze = linspace(dataStruct.dombnd(1,3), ...
                  dataStruct.dombnd(2,3), ...
                  nblockz*nzb+1);
    xc = 0.5*(xe(1:end-1)+xe(2:end));
    yc = 0.5*(ye(1:end-1)+ye(2:end));
    zc = 0.5*(ze(1:end-1)+ze(2:end));
    if terdim
      [dataStruct.uniform.mesh{1} dataStruct.uniform.mesh{2} dataStruct.uniform.mesh{3}] = ndgrid(xc,yc,zc);
    else
      [dataStruct.uniform.mesh{1} dataStruct.uniform.mesh{2}] = ndgrid(xc,yc);
    end

    % compute block interfaces in each direction
    block_intrfc_x = linspace(dataStruct.dombnd(1,1), dataStruct.dombnd(2,1), nblockx+1);
    block_intrfc_y = linspace(dataStruct.dombnd(1,2), dataStruct.dombnd(2,2), nblocky+1);
    block_intrfc_z = linspace(dataStruct.dombnd(1,3), dataStruct.dombnd(2,3), nblockz+1);

    % initialize master arrays
    if precision == 'double'
      for v = 1:dataStruct.nvars
        dataStruct.uniform.(dataStruct.vars{v}) = zeros(nxb*nblockx,nyb*nblocky,nzb*nblockz);
      end
    else
      for v = 1:dataStruct.nvars
        dataStruct.uniform.(dataStruct.vars{v}) = zeros(nxb*nblockx,nyb*nblocky,nzb*nblockz,'single');
      end
    end

    % loop through blocks in data structure
    for blk = 1:dataStruct.nonuniform.nblocks

      % lower ends of the block in each dimension
      xlo = dataStruct.nonuniform.meshbnd(1,1,blk);
      ylo = dataStruct.nonuniform.meshbnd(1,2,blk);
      zlo = dataStruct.nonuniform.meshbnd(1,3,blk);

      % block position in each direction
      block_pos_x = find(block_intrfc_x == xlo);
      block_pos_y = find(block_intrfc_y == ylo);
      block_pos_z = find(block_intrfc_z == zlo);

      % determine starting indices in the master array
      ilo = nxb * (block_pos_x -1) + 1;
      jlo = nyb * (block_pos_y -1) + 1;
      klo = nzb * (block_pos_z -1) + 1;
      ihi = ilo + nxb - 1;
      jhi = jlo + nyb - 1;
      khi = klo + nzb - 1;
%      disp(nxb)
%      disp(nyb)
%      disp(nzb)
%      disp([klo khi])

      % copy data from old to new arrays
      for v = 1:dataStruct.nvars

        % copy
        dataStruct.uniform.(dataStruct.vars{v})(ilo:ihi,jlo:jhi,klo:khi) = dataStruct.nonuniform.(dataStruct.vars{v})(:,:,:,blk);

      end


      % copy mesh data
      %dataStruct.uniform.mesh(ilo:ihi,jlo:jhi,klo:khi+1)

    end

  end

  % delete nonuniform data
  if ~saveNonUni
    for v = 1:dataStruct.nvars
      dataStruct.nonuniform.(dataStruct.vars{v}) = [];
    end
  end

end
