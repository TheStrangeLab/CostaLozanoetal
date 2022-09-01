function [sig,z] = corrDiffTest(r1,n1,r2,n2,alpha)
%CORRDIFFTEST - Test for difference between correlations.
%
% Perform a Fisher (1921) test for a difference between two
% independent correlations (Taken from Howell, David
% C. "Statistical methods for psychology", 5th Edition,
% pp. 277--278). 
%
% If you pass in vectors of equal length, this function will work
% pair-wise to return a vector of significance and Z scores.  The
% significance is determined by seeing if the Z scores is less than
% the absolute value of Z(alpha/2), which makes this a two-tailed
% test.
%
%
% FUNCTION:
%   [sig,z] = corrDiffTest(r1,n1,r2,n2,alpha)
%
% INPUT VARS:
%   r1 - First correlation.
%   n1 - Sample size going into first correlation.
%   r2 - Second correlation.
%   n2 - Sample size of second correlation.
%   alpha - Significance thresh used (defaults to a two-tailed .05,
%           which becomes norminv(alpha/2) in the function).
%
% OUTPUT VARS:
%   sig = 1 or 0 if significantly different or not.
%   z = Z-score for each comparison.
%
%

% process input var
if ~exist('alpha','var')
  alpha = .05;
end

% calculate the z-scores
z = (rprime(r1) - rprime(r2)) ./ sqrt((1./(n1-3))+(1./(n2-3)));

% determine if sig
sig = abs(z) > abs(norminv(alpha/2));



% Calculate an adjusted r
function rp = rprime(r)
rp = .5*log(abs((1+r)./(1-r)));