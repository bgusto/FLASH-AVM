function values = GrabLog(logFile, searchStrings)
%
% GrabLog Get timing data from FLASH log file
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

  % number of strings
  ns = length(searchStrings);

  % open file, read first line
  fid = fopen(logFile);
  tline = fgetl(fid);
  lineCounter = 1;
  while ischar(tline)

    % search for string(s)
    for i = 1:ns
      if length(strfind(tline, searchStrings(i))) > 0
        Key = searchStrings(i);
        Index = strfind(tline, Key);
        values(i) = sscanf(tline(Index(1) + length(Key):end), '%g', 1); 
      end
    end

    % Read next line
    tline = fgetl(fid);
    lineCounter = lineCounter + 1;

  end

  % close
  fclose(fid);

end
