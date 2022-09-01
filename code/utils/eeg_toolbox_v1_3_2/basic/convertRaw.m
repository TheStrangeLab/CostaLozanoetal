function Vout = convertRaw(dat,ampinfo,ampfact)
%CONVERTRAW - Convert raw data to voltage values
%
% FUNCTION:
%   out = convertRaw(dat,ampinfo,ampfact)
%
% INPUT ARGS:
%   dat = rawEEG;              % vector or matrix of raw data
%   ampinfo = [-2048 2048; -5 5];  % Conversion info from raw
%                              %   to voltage (default)
%   ampfact = 10000;           % Amplification factor to
%                              %   correct for (default)
%
% OUTPUT ARGS:
%   Vout - voltage out
%
%


% perform amp conversions
arange = abs(diff(ampinfo(1,:)));
drange = abs(diff(ampinfo(2,:)));

% see if centered around zero for faster operation
if diff(abs(ampinfo(1,:))) == 0 & diff(abs(ampinfo(2,:)))==0
  % are both centered at zero, so get converstion factor
  vfact = (drange*10^6/ampfact)/arange;
  Vout = dat*vfact;
else
  % not centered, so must apply offset
  % adjust for min value in raw analog range
  Vout = dat - min(ampinfo(1,:));
  
  % convert to voltage
  Vout = ((Vout .* drange) ./ arange) + min(ampinfo(2,:));
  
  % divide amp multipler and convert to uV
  Vout = (Vout * 10^6) ./ ampfact;
end

