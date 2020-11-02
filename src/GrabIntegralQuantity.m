function [iq t] = GrabIntegralQuantity(filenm, iqnm, indx)

%-------------------------------------------------------------------------------%
% Info:
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
