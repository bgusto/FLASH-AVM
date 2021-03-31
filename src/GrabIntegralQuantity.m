function [iq t] = GrabIntegralQuantity(filenm, iqnm, indx)
%
% GrabIntegralQuantity Grab data from the integral quantities file
%
%-------------------------------------------------------------------------------%
% Info: This function grabs data from the integral quantities file, written by
%   the IO_writeIntegralQuantities routine of FLASH. The function can take a
%   string argument to match the variable of interest, or the position of the
%   variable in the file can be explicitly specified with the optional argument
%   'indx'. Note that time is always output.
%
% Inputs:
%   filenm - the integral quantities filename
%   iqnm   - the integral quantity name to grab
%   indx   - the index, which can override auto search
%
% Outputs:
%   iq - the integral quantity 
%   t  - time (always output)
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

    % load the data
    [iqraw delim nhead] = importdata(filenm);
    
    % switch if no header detected...
    if nhead == 0

        % check if indx present
        if nargin < 3
            error('If file header not present, user must provide index of desired variable.');
        end

        % get the desired quantitiy + time
        t = iqraw(:,1);
        iq = iqraw(:,indx);

    else

        % provided indx can override search
        if nargin < 3

            % determine desired index by matching string
            for i = 1:length(iqraw.colheaders)
              if contains(iqraw.colheaders{i},iqnm)
                indx = i;
              end
            end

        end

        % get the desired quantitiy + time
        t = iqraw.data(:,1);
        iq = iqraw.data(:,indx);

    end


end
