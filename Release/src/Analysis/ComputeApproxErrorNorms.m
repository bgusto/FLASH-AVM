function dataStruct = ComputeApproxErrorNorms(dataStruct, refDataStruct, var)
%
% ComputeErrorNorms Compute the solution error given the true expression
%
%-------------------------------------------------------------------------------%
% Info:
%
% Inputs:
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

  % initialize norm summations
  l1 = 0.;
  l2 = 0.;
  linf = 0.;

  % loop through the blocks, compute error
  for blk = 1:dataStruct.nonuniform.nblocks

    if dataStruct.ndim == 1

      % evaluate expression at each mesh cell and compute error
      for i = 1:(length(dataStruct.nonuniform.mesh(:,1,blk))-1)

        % compute the bounds of cell i on the adaptive mesh
        xlo = dataStruct.nonuniform.mesh(i,1,blk);
        xhi = dataStruct.nonuniform.mesh(i+1,1,blk);

        % compute cell centers of reference mesh
        cc = 0.5 * (refDataStruct.nonuniform.mesh(1:end-1,:,:) ...
              + refDataStruct.nonuniform.mesh(2:end,:,:));

        % search the reference mesh for cells contained in this bound
        [tcells tblks] = find(cc(:,1,:) > xlo & cc(:,1,:) < xhi);

        % compute the volume-weighted average based on the reference data
        weights = refDataStruct.nonuniform.meshres(1,tblks) / dataStruct.nonuniform.meshres(1,blk); 
        cellavg = 0.;
        for p = 1:length(tcells)
          cellavg = cellavg + weights(p) * refDataStruct.nonuniform.(var)(tcells(p),1,tblks(p));
        end

        % compute residual
        resid = abs(cellavg - dataStruct.nonuniform.(var)(i,1,1,blk));

        % l1 and l2 norms
        l1 = l1 + resid * dataStruct.nonuniform.meshres(1,blk);
        l2 = l2 + resid^2 * dataStruct.nonuniform.meshres(1,blk);
        linf = max(resid, linf);

      end

    elseif dataStruct.ndim == 2

      % evaluate expression at each mesh cell and compute error
      for i = 1:(length(dataStruct.nonuniform.mesh(:,1,blk))-1)
        for j = 1:(length(dataStruct.nonuniform.mesh(1,:,blk))-1)

          % compute the bounds of cell i,j on the adaptive mesh
          xlo = dataStruct.nonuniform.mesh(i,1,blk);
          xhi = dataStruct.nonuniform.mesh(i+1,1,blk);
          ylo = dataStruct.nonuniform.mesh(j,2,blk);
          yhi = dataStruct.nonuniform.mesh(j+1,2,blk);

          % initialize reference cellavg at this target cell
          cellavg = 0.;

          % loop through reference data blocks
          for b0 = 1:refDataStruct.nonuniform.nblocks

            % determine if target cells is within current block's domain
            if xlo >= refDataStruct.nonuniform.mesh(1,1,b0) & ...
                xhi <= refDataStruct.nonuniform.mesh(end,1,b0) & ...
                ylo >= refDataStruct.nonuniform.mesh(1,2,b0) & ...
                yhi <= refDataStruct.nonuniform.mesh(end,2,b0)

              % determine weighting
              weight = (refDataStruct.nonuniform.meshres(1,b0) / dataStruct.nonuniform.meshres(1,blk))^2; 

              % loop through cells in this block
              for i0 = 1:length(refDataStruct.nonuniform.mesh(:,1,1))-1
                for j0 = 1:length(refDataStruct.nonuniform.mesh(:,1,1))-1

                  % determine if current cell center is within the bounds (a subcell of target)
                  ccx = 0.5*(refDataStruct.nonuniform.mesh(i0,1,b0) + refDataStruct.nonuniform.mesh(i0+1,1,b0));
                  ccy = 0.5*(refDataStruct.nonuniform.mesh(j0,2,b0) + refDataStruct.nonuniform.mesh(j0+1,2,b0));
                  if ccx > xlo & ccx < xhi & ccy > ylo & ccy < yhi

                    % compute the volume-weighted average based on the reference data
                    cellavg = cellavg + weight * refDataStruct.nonuniform.(var)(i0,j0,1,b0);

                  end

                end
              end

            end

          end

          % compute residual
          resid = abs(cellavg - dataStruct.nonuniform.(var)(i,j,1,blk));

          % l1 and l2 norms
          l1 = l1 + resid * dataStruct.nonuniform.meshres(1,blk);
          l2 = l2 + resid^2 * dataStruct.nonuniform.meshres(1,blk);
          linf = max(resid, linf);
 
        end
      end

    end

  end

  % take square root of the l2 summation
  l2 = sqrt(l2);

  % add to dataStruct
  dataStruct.nonuniform.(strcat(var, '_l1error')) = l1;
  dataStruct.nonuniform.(strcat(var, '_l2error')) = l2;
  dataStruct.nonuniform.(strcat(var, '_linferror')) = linf;

end
