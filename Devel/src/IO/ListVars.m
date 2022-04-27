function Vars = ListVars(filename)
%
% ListVars Save all variables in the file in a list
%
%-------------------------------------------------------------------------------%
% Info: This function grabs data, and any specified meta-data, from the given
%   hdf5 file. Meta-data might include the node type, bounding box, etc.
%
% Inputs:
%   filename -  the hdf5 filename
%
% Outputs:
%   Vars - a list of variable names in the file
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

  % read data
  hinfo = hdf5info(filename);

  % save all variable names
  for i = 1:length(hinfo.GroupHierarchy.Datasets)
    Vars{i} = erase(hinfo.GroupHierarchy.Datasets(i).Name, '/');
  end

end
