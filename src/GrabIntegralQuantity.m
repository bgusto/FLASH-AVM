function [iq t] = GrabIntegralQuantity(filenm, iqnm)

%-------------------------------------------------------------------------------%
% Info:
%
% Inputs:
%   filenm -  the integral quantities filename
%   ignm -   the integral quantity name to grab
%
% Outputs:
%   iq - the integral quantity 
%   t  - time
%
%-------------------------------------------------------------------------------%

% first parse the header (at least three spaces for delimter)
fid = fopen(filenm);
header = 'header';
nheaders = 1;
while contains(header,'header')
    header = strsplit(fgetl(fid), '   ');
    nheaders = nheaders + 1;
end
fclose(fid);

% determine desired index by matching string
for i = 1:length(header)
  if contains(header{i},iqnm)
    indx = i;
  end
end

% load the data
iqraw = importdata(filenm,' ',nheaders);

% get the desired quantitiy + time
t = iqraw.data(:,1);
iq = iqraw.data(:,indx);
