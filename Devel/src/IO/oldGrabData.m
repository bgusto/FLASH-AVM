function dataStruct = GrabData(filenm, vars, intrp, reflvl)
%
% GrabData Grab data from the specified FLASH checkpoint or plot file.
%
%-------------------------------------------------------------------------------%
% Info:
%   This function grabs data from FLASH hdf5 files for the variable of interest.
%   The function takes block-partitioned data and converts it into a single
%   array. For multi-level adaptive mesh refinement data, the user may choose to
%   prolong all of the data to the finest level using the 'intrp' argument.
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

  % determine dimensionality of dataset
  for blk = 1:nblocks
    if length(tmpData{4}(1,:,1,blk)) > 1, secdim = true;
    if length(tmpData{4}(1,1,:,blk)) > 1, terdim = true;
  end

  % inform data structure of dimensionality
  dataStruct.ndim = 1 + secdim + (secdim & terdim);

  %%%% initialize the global domain boundaries
  %%%xmin = realmax; xmax = realmin;
  %%%if secdim, ymin = realmax; ymax = realmin;
  %%%if terdim, zmin = realmax; zmax = realmin;

  %%%% initialize max level
  %%%maxlvl = 0;

  %%%% initialize minimum mesh spacing (because we want to know what spacing is on finest level)
  %%%dxmin = realmax;
  %%%if secdim, dymin = realmax;
  %%%if terdim, dzmin = realmax;

  % loop through blocks
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
      if secdim, nyb = length(tmpData{4}(1,:,1,blk));
      if terdim, nzb = length(tmpData{4}(1,1,:,blk));

      %%%% determine the bounding box of the global block
      %%%xmin = min(xmin, xlo);
      %%%xmax = max(xmax, xhi);
      %%%if secdim
      %%%  ymin = min(ymin, ylo);
      %%%  ymax = max(ymax, yhi);
      %%%end
      %%%if terdim
      %%%  zmin = min(zmin, zlo);
      %%%  zmax = max(zmax, zhi);
      %%%end

      % cell width
      dx = (xhi - xlo) / nxb;
      if secdim, dy = (yhi - ylo) / nyb;
      if terdim, dz = (zhi - zlo) / nzb;

      % save the solution data in dataStruct
      for i = 1:length(vars)
        dataStruct.nonuniform.(vars{i})(:,:,:,cnt) = tmpData{nmeta+i}(:,:,:,blk);
      end

      % set block's mesh in dataStruct
      dataStruct.nonuniform.mesh(:,cnt,1) = linspace(xlo, xhi, nxb);
      if secdim, dataStruct.nonuniform.mesh(:,cnt,2) = linspace(ylo, yhi, nyb);
      if secdim, dataStruct.nonuniform.mesh(:,cnt,3) = linspace(zlo, zhi, nzb);

      % set mesh block boundaries in dataStruct also
      dataStruct.nonuniform.meshbnd(1,cnt,1) = xlo;
      dataStruct.nonuniform.meshbnd(2,cnt,1) = xhi;
      if secdim
        dataStruct.nonuniform.meshbnd(1,cnt,2) = ylo;
        dataStruct.nonuniform.meshbnd(2,cnt,2) = yhi;
      end
      if terdim
        dataStruct.nonuniform.meshbnd(1,cnt,3) = zlo;
        dataStruct.nonuniform.meshbnd(2,cnt,3) = zhi;
      end

      % add local mesh resolution to dataStruct
      dataStruct.nonuniform.meshres(cnt,1) = dx;
      if secdim, dataStruct.nonuniform.meshres(cnt,2) = dy;
      if terdim, dataStruct.nonuniform.meshres(cnt,3) = dz;

      % add amr refinement level to dataStruct
      dataStruct.nonuniform.lrefine(cnt) = lrefine(blk);

      %%%% keep track of global mesh sizes
      %%%dxmin = min(dxmin, dx);
      %%%if secdim, dymin = min(dymin, dy);
      %%%if terdim, dzmin = min(dzmin, dz);

      %%%% remember which level dxmin found on
      %%%maxlvl = max(maxlvl, lrefine(blk));

    end

    % increment counter
    cnt = cnt + 1;

  end

  % generate output data (prolong or not)
  if prolong

    % determine max level automatically if prolonging and no desired refinement level given
    if nargin < 4
      reflvl = max(lrefine);
    end

    % determine number of base blocks in each direction, create array of block centers
    nbasex = 0;
    if secdim, nbasey = 0;
    if terdim, nbasez = 0;
    blkxcenters = [];
    if secdim, blkycenters = [];
    if terdim, blkzcenters = [];
    for blk = 1:(cnt-1)
      if lrefine(blk) == 1

        % add block center to array
        blkxcenters = [blkxcenters 0.5*(tmpData{2}(1,1,blk)+tmpData{2}(2,1,blk))];
        if secdim, blkycenters = [blkycenters 0.5*(tmpData{2}(1,2,blk)+tmpData{2}(2,2,blk))];
        if terdim, blkzcenters = [blkzcenters 0.5*(tmpData{2}(1,3,blk)+tmpData{2}(2,3,blk))];

      end
    end

    % get unique block centers
    blkxcenters = sort(unique(blkxcenters));
    if secdim, blkycenters = sort(unique(blkycenters));
    if terdim, blkzcenters = sort(unique(blkzcenters));

    % compute number of base blocks based on base block centers
    nbasex = length(blkxcenters);
    if secdim, nbasey = length(blkycenters);
    if terdim, nbasez = length(blkzcenters);

    %%%% create the global domain, compute number of cells
    %%%if reflvl == maxlvl

    %%%  % use the smallest step size in the mesh to build the uniform mesh
    %%%  dataStruct.uniform.mesh(:,1) = xmin:dxmin:xmax;
    %%%  if secdim, dataStruct.uniform.mesh(:,2) = ymin:dymin:ymax;
    %%%  if terdim, dataStruct.uniform.mesh(:,3) = zmin:dzmin:zmax;

    %%%  % number of cells on global mesh
    %%%  nx = length(dataStruct.uniform.mesh(:,1))-1;
    %%%  if secdim, ny = length(dataStruct.uniform.mesh(:,2))-1;
    %%%  if terdim, nz = length(dataStruct.uniform.mesh(:,3))-1;

    %%%else

    % compute number of cells on finest level to cover domain
    nxglob = nbasex * nxb * 2^(reflvl-1);
    if secdim, nyglob = nbasey * nyb * 2^(reflvl-1);
    if terdim, nzglob = nbasez * nzb * 2^(reflvl-1);

    % compute mesh based on number of cells
    dataStruct.uniform.mesh(:,1) = linspace(xmin, xmax, nxglob);
    if secdim, dataStruct.uniform.mesh(:,2) = linspace(ymin, ymax, nyglob);
    if terdim, dataStruct.uniform.mesh(:,3) = linspace(zmin, zmax, nzglob);

    %%%end

    % map the data on leaf blocks to the global data array
    for blk = 1:(cnt-1)

        % ends of block (1=lo,2=hi)
        xlo = tmpData{2}(1,1,blk);
        xhi = tmpData{2}(2,1,blk);

        % get data
        blkdata = tmpData{4}(:,1,1,blk);

        % get amr level
        lvl = lrefine(blk);

        % prolong the data to finest level
        for l = lvl:reflvl-1

          % get size of current block
          nxblk = length(blkdata);

          % create new data for block at higher level
          blkdata2 = zeros(1,2*nxblk);

          % determine how to prolong data
          if strcmp(intrp,'copy')

            % copy parent cells to children
            for i = 1:nxblk
              blkdata2(2*i-1) = blkdata(i);
              blkdata2(2*i) = blkdata(i);
            end

          else

            error('Invalid prolongation option: intrp. Valid options are ... copy')

          end

          % copy variable for next iteration
          blkdata = blkdata2;

        end

        % get size of this newly created block
        nxblk = length(blkdata);

        % determine the starting (i,j) indices in global array
        iglb = int32((xlo - xmin) / dxmin + 1);

        % copy data
        for i = 0:nxblk-1
          vardata(iglb+i) = blkdata(i+1);
        end

      end

    end

  end

end
