function [vdata mdata x y z] = GrabData3D(filenm, varnm, intrp, reflvl)
%
% GrabData3D Grab three-dimensional data from the specified input file.
%
%-------------------------------------------------------------------------------%
% Info: This function grabs data from FLASH hdf5 files for the variable of
%   interest.  The function returns a an array the size of the numbre of blocks,
%   with each element containing that block's structured data (specified by
%   'varnm'). Assumes data is cell-centered (no face-vars). Current prolongation
%   options are 'copy' and 'avginterp.' The former simply copies the parent cell
%   value to the children, while the latter uses a third-order average
%   interpolating polynomial which is biased on the block boundaries.
%
% Inputs:
%   filenm -  the hdf5 filename
%   varnm -   the variable name to be plotted from hdf5 files (a string)
%   intrp  - optional interpolation choice, if prolonging data
%
% Outputs:
%   vdata - the output data matrix
%   mdata - bounding box for each leaf block
%   x - meshgrid, x coordinates
%   y - meshgrid, y coordinates
%   z - meshgrid, z coordinates
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

  % get 'varnm' variable data from current hdf5 data
  Data = GrabHDF5(filenm,{'node type', 'bounding box', 'refine level', varnm});

  % amr metadata 
  nodetyp = Data{1};
  nblocks = length(Data{1});
  lrefine = Data{3};
  nxb = length(Data{4}(:,1,1,1));
  nyb = length(Data{4}(1,:,1,1));
  nzb = length(Data{4}(1,1,:,1));

  % initialize block bounding coordinates array
  mdata = zeros(6,length(find(nodetyp==1)));

  % initialize xmin, xmax, ymin, ymax + dxmin + dymin
  xmin = realmax; xmax = -realmax;
  ymin = realmax; ymax = -realmax;
  zmin = realmax; zmax = -realmax;

  % initialize to large number
  dxmin = realmax;
  dymin = realmax;
  dzmin = realmax;

  % determine max level if prolonging
  if nargin < 4
    reflvl = max(lrefine);
  end

  % determine number of base blocks
  nbase = 0;
  for blk = 1:nblocks
    if lrefine(blk) == 1
      nbase = nbase + 1;
    end
  end

  % initial sweep through data to get xmin, xmax, ymin, ymax + dx, dy
  cnt = 1;
  for blk = 1:nblocks

    % check if leaf block
    if nodetyp(blk) == 1

      % ends of block (1=lo,2=hi)
      xlo = Data{2}(1,1,blk);
      xhi = Data{2}(2,1,blk);
      ylo = Data{2}(1,2,blk);
      yhi = Data{2}(2,2,blk); 
      zlo = Data{2}(1,3,blk);
      zhi = Data{2}(2,3,blk); 

      % determine the bounding box of the global block
      xmin = min(xmin,xlo);
      xmax = max(xmax,xhi);
      ymin = min(ymin,ylo);
      ymax = max(ymax,yhi);
      zmin = min(zmin,zlo);
      zmax = max(zmax,zhi);

      % cell width
      dxmin = min(dxmin, (xhi - xlo) / nxb);
      dymin = min(dymin, (yhi - ylo) / nyb);
      dzmin = min(dzmin, (zhi - zlo) / nzb);

      % save in mesh data
      mdata(1,cnt) = xlo;
      mdata(2,cnt) = xhi;
      mdata(3,cnt) = ylo;
      mdata(4,cnt) = yhi;
      mdata(5,cnt) = zlo;
      mdata(6,cnt) = zhi;

      % increment counter
      cnt = cnt + 1;

    end

  end

  if reflvl == max(lrefine)

    % create the global domain
    x = xmin:dxmin:xmax;
    y = ymin:dymin:ymax;
    z = zmin:dzmin:zmax;

    % number of cells
    nx = length(x)-1;
    ny = length(y)-1;
    nz = length(z)-1;

  else

    % compute number of cells on finest level to cover domain
    nx = sqrt(nbase) * nxb * 2^(reflvl-1);
    ny = sqrt(nbase) * nyb * 2^(reflvl-1);
    nz = sqrt(nbase) * nzb * 2^(reflvl-1);

  end

  % create meshgrid
  [x, y, z] = meshgrid(linspace(xmin,xmax,nx), ...
                        linspace(ymin,ymax,ny), ...
                        linspace(zmin,zmax,nz));

  % initialize the global block
  vdata = zeros(nx,ny,nz);

  % now map the data on blocks to the global data array
  for blk = 1:nblocks

    % check if leaf block
    if nodetyp(blk) == 1

      % ends of block (1=lo,2=hi)
      xlo = Data{2}(1,1,blk);
      xhi = Data{2}(2,1,blk);
      ylo = Data{2}(1,2,blk);
      yhi = Data{2}(2,2,blk); 
      zlo = Data{2}(1,3,blk);
      zhi = Data{2}(2,3,blk); 

      % get data
      blkdata = Data{4}(:,:,:,blk);

      % check amr level
      lvl = lrefine(blk);

      % prolong the data to finest level
      for l = lvl:max(lrefine)-1

        % get size of current block
        [nxblk nyblk nzblk] = size(blkdata);

        % create new data for block at higher level
        blkdata2 = zeros(2*nxblk,2*nyblk,2*nzblk);

        if strcmp(intrp,'copy')

            % loop through cells in block
            for k = 1:nzblk
              for j = 1:nyblk
                for i = 1:nxblk

                  % copy parent cells to children
                  blkdata2(2*i-1,2*j-1,2*k-1) = blkdata(i,j);
                  blkdata2(2*i-1,2*j,2*k-1) = blkdata(i,j);
                  blkdata2(2*i,2*j-1,2*k-1) = blkdata(i,j);
                  blkdata2(2*i,2*j,2*k-1) = blkdata(i,j);
                  blkdata2(2*i-1,2*j-1,2*k) = blkdata(i,j);
                  blkdata2(2*i-1,2*j,2*k) = blkdata(i,j);
                  blkdata2(2*i,2*j-1,2*k) = blkdata(i,j);
                  blkdata2(2*i,2*j,2*k) = blkdata(i,j);

                end
              end
            end

        else

            error('Invalid prolongation option: intrp. Valid options are: copy')

        end

        % copy variable
        blkdata = blkdata2;

      end

      % get size of this newly created block
      [nxblk nyblk nzblk] = size(blkdata);

      % determine the starting (i,j) indices in global array
      iglb = int32((xlo - xmin) / dxmin + 1);
      jglb = int32((ylo - ymin) / dymin + 1);
      kglb = int32((zlo - zmin) / dzmin + 1);

      % copy data
      for k = 0:nzblk-1
        for j = 0:nyblk-1
          for i = 0:nxblk-1
            vdata(iglb+i,jglb+j,kglb+k) = blkdata(i+1,j+1,k+1);
          end
        end
      end

    end

  end

end
