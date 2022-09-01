function [pthresh_pos,pthresh_neg] = get_pthresh(t1_desired,p_boot,combineBins)
%GET_PTHRESH = Gets threshold for desired Type-I from bootstrap.
%
% Calculates the p thresholds from a bootstrap analysis.
% 
% p_boot is made up as:
%
%   p_boot(iterations,bin1,bin2,bin3)  
%
% FUNCTION:
%   [pthresh_pos,pthresh_neg] = get_pthresh(t1_desired,p_boot,combineBins)
%
% INPUT ARGUMENTS:
%   t1_desired = .01; % Percent of leads by chance
%   p_boot = p_boot(:,goodleads,:); % from the bootstrap
%   combineBins = 0;  % 1 to combine all bins to a big distribution
%                     % defaults to 0 for backwards compat.
%
% OUTPUT ARGS:
%   pthresh_pos = pthresh(bin1,bin2,bin3); % the positive p_thresh
%                                     % for bin1 & bin2
%
%   pthresh_neg = pthresh(bin1,bin2,bin3); % the positive p_thresh
%                                     % for bin1 & bin2
%
%   If combineBins==1, then pthresh will be single values that you
%   should apply to all bins of your data.
%

% 2004/07/03 - PBS - Removed squeeze from return vals.

if ~exist('combineBins','var')
  combineBins = 0;
end

% calculate the size of the distribution
if combineBins
  % combine all dimensions
  distSize = prod(size(p_boot));
else
  % only use iterations
  distSize = size(p_boot,1);
end

% calc the thresh
t1 = fix(distSize*t1_desired)+1;

if combineBins
  % reshape and sort the entire thing
  p_boot = sort(reshape(p_boot,distSize,1));
  
  % get the thresh
  pthresh_pos = p_boot(t1);
  pthresh_neg = 1-(p_boot(distSize-t1+1));
else
  % loop over each dimension and sort the iterations
  for b1 = 1:size(p_boot,2)
    for b2 = 1:size(p_boot,3)
      for b3 = 1:size(p_boot,4)
        
        % sort the iterations
        p_boot(:,b1,b2,b3) = sort(squeeze(p_boot(:,b1,b2,b3)));
        
      end
    end
  end
  
  % get the thresh
  pthresh_pos = shiftdim(p_boot(t1,:,:,:));
  pthresh_neg = shiftdim(1-p_boot(size(p_boot,1)-t1+1,:,:,:));
  %sizepb = size(p_boot);
  %pthresh_pos = reshape(double(p_boot(t1,:,:,:)),sizepb(2:end));
  %pthresh_neg = reshape(1 - double(p_boot(size(p_boot,1)-t1+1,:,:,:)),sizepb(2:end));
end





