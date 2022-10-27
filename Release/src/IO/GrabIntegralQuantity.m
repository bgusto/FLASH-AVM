function data = GrabIntegralQuantity(filenm, iqnm)
%
% GrabIntegralQuantity Grab data from the integral quantities file
%
%-------------------------------------------------------------------------------%
% Info: This function grabs data from the integral quantities file, written by
%   the IO_writeIntegralQuantities routine of FLASH. The function retrieves all
%   IO_writeIntegralQuantities data and stores in a structure called data.
%
% Inputs:
%   filenm - the integral quantities filename
%
% Outputs:
%   data - the struct containing headers and data
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

  % load the data with importdata
  raw = importdata(filenm);

  % get the quantity names (eliminates single whitespace between words (ex: 'filling factor' -> 'filling_factor')
  quantities = split(regexprep(strtrim(raw.textdata), '(\w) (\w)', '$1_$2'));

  % store number of quantities found
  nq = length(quantities);

  % eliminate any incompatible characters here (add your own? sometimes IO_writeIntegralQuantities is user modified)
  for i = 1:nq
    quantities{i} = erase(quantities{i}, {'#', '-'});
  end

  % organize raw data by quantity
  if nargin > 1 

    for j = 1:length(iqnm)
      for i = 1:nq
        if strcmp(quantities{i},iqnm{j})
          data.(quantities{i}) = raw.data(:,i);
        end
      end
    end

  else

    for i = 1:nq
      data.(quantities{i}) = raw.data(:,i);
    end

  end

end
