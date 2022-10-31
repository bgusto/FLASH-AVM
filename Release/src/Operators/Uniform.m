function dataStruct = Uniform(dataStruct, refLevel, interpType, saveNonUni, precision)
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
  if nargin < 2 | refLevel == 0
    refLevel = dataStruct.lrefine_max;
  end
  if nargin < 3
    interpType = 'copy';
  end
  if nargin < 4
    saveNonUni = true;
  end
  if nargin < 5
    precision = 'double';
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

  % difference in levels from finest to desired 
  lref_max = max(dataStruct.nonuniform.lrefine);
  dl = double(refLevel) - double(lref_max);

  % abort if nxb*2^(dl) < 1
  if nxb*2^(dl) < 1 | nyb*2^(dl) < 1
    error('Cannot coarsen the block data to the desired level. Raise desired refinement level.')
  end

  % get the cell resolution at the desired level
  dx_desired = 2^(-dl) * min(dataStruct.nonuniform.meshres(1,:));
  dy_desired = 2^(-dl) * min(dataStruct.nonuniform.meshres(2,:));
  if terdim
    dz_desired = 2^(-dl) * min(dataStruct.nonuniform.meshres(3,:));
  end

  % create the uniform mesh
  xe = xlo:dx_desired:xhi;
  ye = ylo:dy_desired:yhi;
  if terdim
    ze = zlo:dz_desired:zhi;
  end
  xc = 0.5*(xe(1:end-1)+xe(2:end));
  nx = length(xc);
  yc = 0.5*(ye(1:end-1)+ye(2:end));
  ny = length(yc);
  if terdim
    zc = 0.5*(ze(1:end-1)+ze(2:end));
    nz = length(zc);
  else
    nz = 1;
  end
  if terdim
    [dataStruct.uniform.mesh{1} dataStruct.uniform.mesh{2} dataStruct.uniform.mesh{3}] = ndgrid(xc,yc,zc);
  else
    [dataStruct.uniform.mesh{1} dataStruct.uniform.mesh{2}] = ndgrid(xc,yc);
  end

  % initialize master arrays
  if strcmp(precision, 'double')
    for v = 1:dataStruct.nvars
      dataStruct.uniform.(dataStruct.vars{v}) = zeros(nx,ny,nz);
    end
  else
    for v = 1:dataStruct.nvars
      dataStruct.uniform.(dataStruct.vars{v}) = zeros(nx,ny,nz,'single');
    end
  end

  % loop through nonuniform structure and prolong if necessary to uniform grid
  for blk = 1:dataStruct.nonuniform.nblocks

    % lower ends of the block in each dimension
    bxlo = dataStruct.nonuniform.meshbnd(1,1,blk);
    bxhi = dataStruct.nonuniform.meshbnd(2,1,blk);
    bylo = dataStruct.nonuniform.meshbnd(1,2,blk);
    byhi = dataStruct.nonuniform.meshbnd(2,2,blk);
    bzlo = dataStruct.nonuniform.meshbnd(1,3,blk);
    bzhi = dataStruct.nonuniform.meshbnd(2,3,blk);

    % find where tmp fits into uniform grid at the finest level
    ilo = find(xe == bxlo);
    ihi = find(xe == bxhi);
    jlo = find(ye == bylo);
    jhi = find(ye == byhi);
    if terdim
      klo = find(ze == bzlo);
      khi = find(ze == bzhi);
    end

    % difference between desired level and current block's level
    dl = refLevel - dataStruct.nonuniform.lrefine(blk);

    % prolong up to finest level (or coarsen)
    if dl > 0

      % loop through vars
      for v = 1:dataStruct.nvars

        % temporary data field
        if terdim
          tmp = zeros((2^dl)*nxb,(2^dl)*nyb,(2^dl)*nzb);
        else
          tmp = zeros((2^dl)*nxb,(2^dl)*nyb,1);
        end

        % copy data or interpolate
        if interpType == 'copy'

          if terdim
            for k = 1:nzb
              for j = 1:nyb
                for i = 1:nxb
                  tmp(i,j,k) = dataStruct.nonuniform.(dataStruct.vars{v})(:,:,:,blk);
                end
              end
            end
          else
            for j = 1:nyb
              for i = 1:nxb
                tmp((2^dl*(i-1)+1):(2^dl*(i-1)+2^dl),(2^dl*(j-1)+1):(2^dl*(j-1)+2^dl)) = dataStruct.nonuniform.(dataStruct.vars{v})(i,j,1,blk) * ones(2^dl,2^dl);
              end
            end
          end

        end

        % copy tmp data to new uniform grid datastructure
        if terdim
          dataStruct.uniform.(dataStruct.vars{v})(ilo:(ilo+2^dl*nxb-1),jlo:(jlo+2^dl*nyb-1),klo:(klo+2^dl*nzb-1)) = tmp;
        else
          dataStruct.uniform.(dataStruct.vars{v})(ilo:(ilo+2^dl*nxb-1),jlo:(jlo+2^dl*nyb-1),1) = tmp;
        end

      end

    else % coarsen

      % loop through vars
      for v = 1:dataStruct.nvars

        % temporary data storage
        if terdim
          tmpFine = dataStruct.nonuniform.(dataStruct.vars{v})(:,:,:,blk);
        else
          tmpFine = dataStruct.nonuniform.(dataStruct.vars{v})(:,:,1,blk);
        end

        % loop until coarsest level
        for l = 1:abs(dl)-1

          % initialize coarse variable
          if terdim
            tmpCoarse = zeros(nxb/2^(-l),nyb/2^(-l),nzb/2^(-l));
          else
            tmpCoarse = zeros(nxb/2^(-l),nyb/2^(-l),1);
          end

          if terdim
            for k = 1:nzb/2^(-l)
              for j = 1:nyb/2^(-l)
                for i = 1:nxb/2^(-l)
                  tmpCoarse(i,j,k) = 0.125*(tmpFine(2*i-1,2*j-1,2*k-1) + ...
                                            tmpFine(2*i-1,2*j-1,2*k) + ...
                                            tmpFine(2*i-1,2*j,2*k-1) + ...
                                            tmpFine(2*i,2*j-1,2*k-1) + ...
                                            tmpFine(2*i-1,2*j,2*k) + ...
                                            tmpFine(2*i,2*j-1,2*k) + ...
                                            tmpFine(2*i,2*j,2*k-1) + ...
                                            tmpFine(2*i,2*j,2*k));
                end
              end
            end
          else
            for j = 1:nyb/2^(-l)
              for i = 1:nxb/2^(-l)
                tmpCoarse(i,j,k) = 0.25*(tmpFine(2*i,2*j-1,1) + ...
                                          tmpFine(2*i-1,2*j,1) + ...
                                          tmpFine(2*i,2*j,1) + ...
                                          tmpFine(2*i-1,2*j-1,1));
              end
            end
          end

          % reset variables
          tmpFine = tmpCoarse;

        end

        % copy tmp data to new uniform grid datastructure
        if terdim
          dataStruct.uniform.(dataStruct.vars{v})(ilo:(ilo+2^dl*nxb-1),jlo:(jlo+2^dl*nyb-1),klo:(klo+2^dl*nzb-1)) = tmpFine;
        else
          dataStruct.uniform.(dataStruct.vars{v})(ilo:(ilo+2^dl*nxb-1),jlo:(jlo+2^dl*nyb-1),1) = tmpFine;
        end

      end

    end

  end

  % delete nonuniform data
  if ~saveNonUni
    for v = 1:dataStruct.nvars
      dataStruct.nonuniform.(dataStruct.vars{v}) = [];
    end
  end

end
