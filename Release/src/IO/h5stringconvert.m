function strs = h5stringconvert(h5strings)
%H5STRINGCONVERT converts strings read from HDF5 files to char arrays.
% STRS = H5STRINGCONVERT(H5STRINGS) converts strings read from an HDF5 file
% by the hdf5read command (H5STRINGS) to a string (char array) or cell
% array of strings if there are multiple, which is returned in STRS. This
% function is necessary because the two different ways that strings can be
% stored in HDF5 files yield different return types from hdf5read.

% Copyright (c) 2013, Freja Nordsiek
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met
%
% 1. Redistributions of source code must retain the above copyright notice, this
% list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% It is either already a string or an hdf5.h5string. If the former, it
% can be returned as is. If the latter, then it must be extracted.

if ischar(h5strings)

    % Its already a string, so simply return it.
    strs = h5strings;

elseif numel(h5strings) == 1 && isa(h5strings,'hdf5.h5string')

    % It is a single hdf5.h5string, so extract it.
    strs = h5strings.Data;

elseif iscell(h5strings) || isa(h5strings,'hdf5.h5string')

    % It is a cell array or an array of hdf5.h5string's. So, we make a
    % cell array for strs and convert entry by entry by a recursive
    % call (cell array) or extraction (hdf5.h5string).

    strs = cell(size(h5strings));

    for ii=1:numel(h5strings)
        if iscell(h5strings)
            strs{ii} = hdf5stringconvert(h5strings{ii});
        else
            strs{ii} = h5strings(ii).Data;
        end
    end

else
    error('Not a string.');
end
