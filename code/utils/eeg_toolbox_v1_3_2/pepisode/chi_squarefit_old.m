function [thresh,R2] = chi_squarefit_old(freqs, pows)
% calc the fit for the pepisode
nbins=1000;

pv=polyfit(log10(freqs),log10(pows'),1); % linear regression
means=10.^(polyval(pv,log10(freqs))); % these are the fitted mean powers

R=corrcoef(polyval(pv,log10(freqs)),log10(pows'));
R2=R(1,2)^2;

% Now, compute the chi2cdf vals.
thresh=[freqs;chi2inv(0:(1/nbins):(1-(1/nbins)),2)'*(means/2)];
% divide the means by two because the mean of the chi-square distribution is equal to the computed mean divided by the degrees of freedom
