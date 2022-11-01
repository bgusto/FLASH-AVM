function dataStruct = ComputeErrorNorms(dataStruct, var, expression, intmethod)
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

        % compute volume integral using fancy quadrature
        %cellavg = integral(expression, ...
        %              dataStruct.nonuniform.mesh(i,1,blk), ...
        %              dataStruct.nonuniform.mesh(i+1,1,blk)); ...
                      %/ dataStruct.nonuniform.meshres(1,blk);

        % choose integration method
        if strcmp(intmethod, 'midpoint')

          cellavg = expression(...
                        0.5*(dataStruct.nonuniform.mesh(i,1,blk) ...
                        + dataStruct.nonuniform.mesh(i+1,1,blk)));

        elseif strcmp(intmethod, 'trapezoid')

          % number of subintervals
          nsubint = 8;

          % subinterval spacing
          dxsub = abs(dataStruct.nonuniform.mesh(i,1,blk) ...
                        - dataStruct.nonuniform.mesh(i+1,1,blk)) / nsubint;

          % submesh
          xsub = linspace(dataStruct.nonuniform.mesh(i,1,blk), ...
                  dataStruct.nonuniform.mesh(i+1,1,blk), nsubint+1);

          % compute integral
          cellavg = 0.5 * dxsub * sum(expression(xsub(1:end-1)) + expression(xsub(2:end))) ...
                     / dataStruct.nonuniform.meshres(1,blk);

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
        for j = 1:(length(dataStruct.nonuniform.mesh(:,2,blk))-1)

          % choose integration method
          if strcmp(intmethod, 'midpoint')

            % compute volume integral average using midpoint rule
            cellavg = expression(...
                          0.5*(dataStruct.nonuniform.mesh(i,1,blk) ...
                          + dataStruct.nonuniform.mesh(i+1,1,blk)), ...
                          0.5*(dataStruct.nonuniform.mesh(j,2,blk) ...
                          + dataStruct.nonuniform.mesh(j+1,2,blk)));

          elseif strcmp(intmethod, 'trapezoid')

            % number of subintervals per side
            nsubint = 8;

            % subinterval spacing
            dxsub = abs(dataStruct.nonuniform.mesh(i,1,blk) ...
                          - dataStruct.nonuniform.mesh(i+1,1,blk)) / nsubint;
            dysub = abs(dataStruct.nonuniform.mesh(j,2,blk) ...
                          - dataStruct.nonuniform.mesh(j+1,2,blk)) / nsubint;

            % submesh
            xsub = linspace(dataStruct.nonuniform.mesh(i,1,blk), ...
                    dataStruct.nonuniform.mesh(i+1,1,blk), nsubint+1);
            ysub = linspace(dataStruct.nonuniform.mesh(j,2,blk), ...
                    dataStruct.nonuniform.mesh(j+1,2,blk), nsubint+1);

            % add edge values to summation
            sm = expression(xsub(1), ysub(1)) ...
                  + expression(xsub(end), ysub(1)) ...
                  + expression(xsub(1), ysub(end)) ...
                  + expression(xsub(end), ysub(end));

            for ii = 2:length(xsub)-1
              sm = sm + 2*expression(xsub(ii), ysub(1));
            end
            for ii = 2:length(xsub)-1
              sm = sm + 2*expression(xsub(ii), ysub(end));
            end
            for jj = 2:length(ysub)-1
              sm = sm + 2*expression(xsub(1), ysub(jj));
            end
            for jj = 2:length(ysub)-1
              sm = sm + 2*expression(xsub(end), ysub(jj));
            end
            for ii = 2:length(xsub)-1
              for jj = 2:length(ysub)-1
                sm = sm + 4*expression(xsub(ii), ysub(jj));
              end
            end
            sm = sm * dxsub * dysub / 4;
            cellavg = sm / (dataStruct.nonuniform.meshres(1,blk) ...
                  * dataStruct.nonuniform.meshres(2,blk));

          elseif strcmp(intmethod, 'default')

            % fancy quadrature
            cellavg = integral2(expression, ...
                          dataStruct.nonuniform.mesh(i,1,blk), ...
                          dataStruct.nonuniform.mesh(i+1,1,blk), ...
                          dataStruct.nonuniform.mesh(j,2,blk), ...
                          dataStruct.nonuniform.mesh(j+1,2,blk)) ...
                          / (dataStruct.nonuniform.meshres(1,blk) ...
                          * dataStruct.nonuniform.meshres(2,blk));

          end

          % compute residual
          resid = abs(cellavg - dataStruct.nonuniform.(var)(i,j,1,blk));
          %dataStruct.nonuniform.(strcat(var, '_resid'))(i,j,blk) = resid;

          % l1 and l2 norms
          l1 = l1 + resid * dataStruct.nonuniform.meshres(1,blk) ...
                * dataStruct.nonuniform.meshres(2,blk);
          l2 = l2 + resid^2 ...
                * dataStruct.nonuniform.meshres(1,blk) ...
                * dataStruct.nonuniform.meshres(2,blk);
          linf = max(resid, linf);

        end
      end

    else

    end

  end

  % take square root of the l2 summation
  l2 = sqrt(l2);

  % add to dataStruct
  dataStruct.nonuniform.(strcat(var, '_l1error')) = l1;
  dataStruct.nonuniform.(strcat(var, '_l2error')) = l2;
  dataStruct.nonuniform.(strcat(var, '_linferror')) = linf;

end
