function dataStruct = GrabData(filenm, vars, precision)
%
% GrabData Grab data from the specified FLASH checkpoint or plot file.
%
%-------------------------------------------------------------------------------%
% Info:
%   This function grabs data from FLASH hdf5 files for the variable of interest.
%   The function assumes a Paramesh tree structure.
%
% Inputs:
%   filenm - the hdf5 filename
%   vars   - cell array of desired output variables (e.g. {'dens', 'pres', 'velx'}, or just {'dens'})
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

  % default arguments
  if nargin < 3
    precision = 'single';
  end

  % construct cell array to pass
  inputs = [{'node type', 'bounding box', 'refine level', 'real scalars', 'integer scalars'} vars];

  % get 'varnm' variable data from current hdf5 data
  tmpData = GrabHDF5(filenm, inputs);

  % number of meta data elements
  nmeta = 5;

  % get simulation time
  for m = 1:length(tmpData{4})
    scalarnm = strtrim(h5stringconvert(tmpData{4}(m).Data{1}));
    if strcmp(erase(scalarnm,{''''}), 'time')
      dataStruct.simtime = tmpData{4}(m).Data{2};
    elseif strcmp(erase(scalarnm,{''''}), 'dt')
      dataStruct.simdt = tmpData{4}(m).Data{2};
    end
  end

  % get nstep
  for m = 1:length(tmpData{5})
    scalarnm = strtrim(h5stringconvert(tmpData{5}(m).Data{1}));
    if strcmp(erase(scalarnm,{''''}), 'nstep')
      dataStruct.nstep = tmpData{5}(m).Data{2};
    end
  end

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
    if length(tmpData{nmeta+1}(1,:,1,blk)) > 1, secdim = true; end
    if length(tmpData{nmeta+1}(1,1,:,blk)) > 1, terdim = true; end
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
      nxb = length(tmpData{nmeta+1}(:,1,1,blk));
      if secdim, nyb = length(tmpData{nmeta+1}(1,:,1,blk)); end
      if terdim, nzb = length(tmpData{nmeta+1}(1,1,:,blk)); end

      % save nxb, nyb, nzb
      data.nonuniform.ncb = [nxb];
      if secdim, data.nonuniform.ncb = [nxb nyb]; end
      if terdim, data.nonuniform.ncb = [nxb nyb nzb]; end

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
      else
        dataStruct.nonuniform.meshbnd(1,2,cnt) = 0;
        dataStruct.nonuniform.meshbnd(2,2,cnt) = 0;
      end
      if terdim
        dataStruct.nonuniform.meshbnd(1,3,cnt) = zlo;
        dataStruct.nonuniform.meshbnd(2,3,cnt) = zhi;
      else
        dataStruct.nonuniform.meshbnd(1,3,cnt) = 0;
        dataStruct.nonuniform.meshbnd(2,3,cnt) = 0;
      end

      % save the solution data in dataStruct
      if precision == 'double'
        for i = 1:dataStruct.nvars
          dataStruct.nonuniform.(vars{i})(:,:,:,cnt) = tmpData{nmeta+i}(:,:,:,blk);
        end
      else
        for i = 1:dataStruct.nvars
          dataStruct.nonuniform.(vars{i})(:,:,:,cnt) = single(tmpData{nmeta+i}(:,:,:,blk));
        end
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
  dataStruct.dombnd(2,1) = xmax;
  if secdim
    dataStruct.dombnd(1,2) = ymin;
    dataStruct.dombnd(2,2) = ymax;
  else
    dataStruct.dombnd(1,2) = 0;
    dataStruct.dombnd(2,2) = 0;
  end
  if terdim
    dataStruct.dombnd(1,3) = zmin;
    dataStruct.dombnd(2,3) = zmax;
  else
    dataStruct.dombnd(1,3) = 0;
    dataStruct.dombnd(2,3) = 0;
  end

end
