function dataStruct = GrabData(filenm, vars)
%
% GrabData Grab data from the specified FLASH checkpoint or plot file.
%
%-------------------------------------------------------------------------------%
% Info:
%   This function grabs data from FLASH hdf5 files for the variable of interest.
%   The function takes block-partitioned data and converts it into a single
%   array. For multi-level adaptive mesh refinement data, the user may choose to
%   prolong all of the data to the finest level using the 'intrp' argument. The
%   function assumes a Paramesh tree structure.
%
% Inputs:
%   filenm - the hdf5 filename
%   vars   - cell array of desired output variables (e.g. {'dens', 'pres', 'velx'}, or just {'dens'})
%   intrp  - optional interpolation choice, if prolonging data
%   reflvl - the desired refinement level, using the coarsest parent block
%            found in the input file as the base level (found automatically if
%            not provided and compatible 'intrp' option specified)
%
% Outputs:
%   dataStruct - structure containing all data fields (solution field, mesh data, etc.)
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

  % determine if we want to prolong AMR data to finest level
  if nargin < 3
    prolong = false;
  else
    prolong = true;
  end

  % construct cell array to pass
  inputs = [{'node type', 'bounding box', 'refine level'} vars];

  % get 'varnm' variable data from current hdf5 data
  tmpData = GrabHDF5(filenm, inputs);

  % number of meta data elements
  nmeta = 3;

  % get some mesh metadata
  nodetyp = tmpData{1};                  % the type of AMR block (leaf or not)
  nblocks = length(tmpData{1});          % the number of AMR blocks in the mesh
  lrefine = tmpData{3};                  % the refinement level of the block ID

  % save the number of variables + name strings
  dataStruct.nvars = length(vars);
  for i = 1:dataStruct.nvars
    vars{i} = erase(vars{i}, ' ');
    dataStruct.vars{i} = vars{i};
  end

  % determine dimensionality of dataset
  secdim = false;
  terdim = false;
  for blk = 1:nblocks
    if length(tmpData{4}(1,:,1,blk)) > 1, secdim = true; end
    if length(tmpData{4}(1,1,:,blk)) > 1, terdim = true; end
  end

  % inform data structure of dimensionality
  dataStruct.ndim = 1 + secdim + (secdim & terdim);

  % initialize the global domain boundaries
  xmin = realmax; xmax = realmin;
  if secdim, ymin = realmax; ymax = realmin; end
  if terdim, zmin = realmax; zmax = realmin; end

  % loop through all blocks
  cnt = 1;
  for blk = 1:nblocks

    % check if leaf block
    if nodetyp(blk) == 1

      % ends of block (1=lo,2=hi)
      xlo = tmpData{2}(1,1,blk);
      xhi = tmpData{2}(2,1,blk);
      ylo = tmpData{2}(1,2,blk);
      yhi = tmpData{2}(2,2,blk);
      zlo = tmpData{2}(1,3,blk);
      zhi = tmpData{2}(2,3,blk);

      % get number of cells in block
      nxb = length(tmpData{4}(:,1,1,blk));
      if secdim, nyb = length(tmpData{4}(1,:,1,blk)); end
      if terdim, nzb = length(tmpData{4}(1,1,:,blk)); end

      % cell width
      dx = (xhi - xlo) / nxb;
      if secdim, dy = (yhi - ylo) / nyb; end
      if terdim, dz = (zhi - zlo) / nzb; end

      % add local mesh resolution to dataStruct
      dataStruct.nonuniform.meshres(1,cnt) = dx;
      if secdim, dataStruct.nonuniform.meshres(2,cnt) = dy; end
      if terdim, dataStruct.nonuniform.meshres(3,cnt) = dz; end

      % add amr refinement level to dataStruct
      dataStruct.nonuniform.lrefine(cnt) = lrefine(blk);

      % set block's mesh in dataStruct
      dataStruct.nonuniform.mesh(1:nxb+1,1,cnt) = linspace(xlo, xhi, nxb+1);
      if secdim, dataStruct.nonuniform.mesh(1:nyb+1,2,cnt) = linspace(ylo, yhi, nyb+1); end
      if terdim, dataStruct.nonuniform.mesh(1:nzb+1,3,cnt) = linspace(zlo, zhi, nzb+1); end

      % set mesh block boundaries in dataStruct also
      dataStruct.nonuniform.meshbnd(1,1,cnt) = xlo;
      dataStruct.nonuniform.meshbnd(2,1,cnt) = xhi;
      if secdim
        dataStruct.nonuniform.meshbnd(1,2,cnt) = ylo;
        dataStruct.nonuniform.meshbnd(2,2,cnt) = yhi;
      end
      if terdim
        dataStruct.nonuniform.meshbnd(1,3,cnt) = zlo;
        dataStruct.nonuniform.meshbnd(2,3,cnt) = zhi;
      end

      % save the solution data in dataStruct
      for i = 1:dataStruct.nvars
        dataStruct.nonuniform.(vars{i})(:,:,:,cnt) = tmpData{nmeta+i}(:,:,:,blk);
      end

      % keep track of global min and max in case prolongation is to be used
      xmin = min(xmin, xlo);
      xmax = max(xmax, xhi);
      if secdim
        ymin = min(ymin, ylo);
        ymax = max(ymax, yhi);
      end
      if terdim
        zmin = min(zmin, zlo);
        zmax = max(zmax, zhi);
      end

      % increment counter
      cnt = cnt + 1;

    end

  end

  % save the number of leaf blocks
  dataStruct.nonuniform.nblocks = cnt - 1;

  % save the global xmin, xmax, ...
  dataStruct.dombnd(1,1) = xmin;
  dataStruct.dombnd(2,1) = xmin;
  if secdim
    dataStruct.dombnd(1,2) = ymin;
    dataStruct.dombnd(2,2) = ymin;
  end
  if terdim
    dataStruct.dombnd(1,3) = zmin;
    dataStruct.dombnd(2,3) = zmin;
  end
  % generate output data (prolong or not)
  %%%if prolong

  %%%  % determine max level automatically if prolonging and no desired refinement level given
  %%%  if nargin < 4
  %%%    reflvl = max(lrefine);
  %%%  end

  %%%  % determine number of base blocks in each direction, create array of block centers
  %%%  nbasex = 0;
  %%%  if secdim, nbasey = 0;
  %%%  if terdim, nbasez = 0;
  %%%  blkxcenters = [];
  %%%  if secdim, blkycenters = [];
  %%%  if terdim, blkzcenters = [];
  %%%  for blk = 1:nblocks
  %%%    if lrefine(blk) == 1

  %%%      % add block center to array
  %%%      blkxcenters = [blkxcenters 0.5*(tmpData{2}(1,1,blk)+tmpData{2}(2,1,blk))];
  %%%      if secdim, blkycenters = [blkycenters 0.5*(tmpData{2}(1,2,blk)+tmpData{2}(2,2,blk))];
  %%%      if terdim, blkzcenters = [blkzcenters 0.5*(tmpData{2}(1,3,blk)+tmpData{2}(2,3,blk))];

  %%%    end
  %%%  end

  %%%  % get unique block centers
  %%%  blkxcenters = sort(unique(blkxcenters));
  %%%  if secdim, blkycenters = sort(unique(blkycenters));
  %%%  if terdim, blkzcenters = sort(unique(blkzcenters));

  %%%  % compute number of base blocks based on base block centers
  %%%  nbasex = length(blkxcenters);
  %%%  if secdim, nbasey = length(blkycenters);
  %%%  if terdim, nbasez = length(blkzcenters);

  %%%  % compute number of cells on finest level to cover domain
  %%%  nxglob = nbasex * nxb * 2^(reflvl-1);
  %%%  if secdim, nyglob = nbasey * nyb * 2^(reflvl-1);
  %%%  if terdim, nzglob = nbasez * nzb * 2^(reflvl-1);

  %%%  % compute mesh based on number of cells
  %%%  dataStruct.uniform.mesh(:,1) = linspace(xmin, xmax, nxglob+1);
  %%%  if secdim, dataStruct.uniform.mesh(:,2) = linspace(ymin, ymax, nyglob+1);
  %%%  if terdim, dataStruct.uniform.mesh(:,3) = linspace(zmin, zmax, nzglob+1);

  %%%  % map the data on leaf blocks to the global data array
  %%%  for blk = 1:(cnt-1)

  %%%    % get number of cells in each direction
  %%%    nxb = length(dataStruct.nonuniform.mesh(:,blk,1));
  %%%    if secdim
  %%%      nyb = length(dataStruct.nonuniform.mesh(:,blk,2));
  %%%    else
  %%%      nyb = 0;
  %%%    end
  %%%    if terdim
  %%%      nzb = length(dataStruct.nonuniform.mesh(:,blk,3));
  %%%    else
  %%%      nzb = 0;
  %%%    end

  %%%    % loop through user-selected variables
  %%%    for v = 1:length(vars)

  %%%      % make temporary copy of solution data
  %%%      blkdata = dataStruct.nonuniform.(vars{v})(:,:,:,blk);

  %%%      % prolong the data to finest level
  %%%      for l = dataStruct.nonuniform.lrefine(blk):reflvl-1

  %%%        % create new temporary storage space for block at higher level
  %%%        blkdata2 = zeros(2*nxb, max(1,2*nyb), max(1,2*nzb));

  %%%        % determine how to prolong data
  %%%        if strcmp(intrp,'copy')

  %%%          switch dataStruct.ndim

  %%%            case 1

  %%%              % copy parent cells to children
  %%%              for i = 1:nxb
  %%%                blkdata2(2*i-1) = blkdata(i,1,1);
  %%%                blkdata2(2*i) = blkdata(i,1,1);
  %%%              end

  %%%            end

  %%%          else

  %%%            error('Invalid prolongation option: intrp. Valid options are ... copy')

  %%%          end

  %%%          % copy variable for next iteration
  %%%          blkdata = blkdata2;

  %%%        end

  %%%        % get size of this newly created block
  %%%        nxblk = length(blkdata);

  %%%        % determine the starting (i,j) indices in global array
  %%%        iglb = int32((xlo - xmin) / dxmin + 1);

  %%%        % copy data
  %%%        for i = 0:nxblk-1
  %%%          vardata(iglb+i) = blkdata(i+1);
  %%%        end

  %%%      end

  %%%  end

  %%%end

end
