function [p, p_z, appoi, apcd_A, apcd_B, surrogate] = ppoi(data1, data2, nperm)
% This is a simple adaptation of Rufin's PhaseOpposition function (see
% VanRullen 2016 Front Neurosci) to allow for vector strength comparisons
% between two experimental conditions that come from phase-to-amplitude
% coupling measures (see Axmacher et al 2010 PNAS)
%
% Diego Lozano-Soldevilla, CTB-UPM, 08-Nov-2019 09:29:50
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright Rufin VanRullen, 2016%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This and associated functions are provided without explicit or implied
%guarantee, and without any form of technical support.
%The function is free to use for any personal, scientific or commercial
%application. When using this function in any published study, please
%cite: "VanRullen, R. (2016). How to evaluate phase differences between
%trial groups in ongoing electrophysiological signals. (submitted)"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Usage: [p_circWW, p_POS, p_zPOS] = P/home/dieloz/Dropbox/Velarde2019Neuroimage/Time-Locked-Index-master/v0haseOpposition(data1, data2, [nperm], [circww_ITCthreshold])
%Computes p-values of phase opposition for a dataset using 2 methods, the
%circular Watson-Williams test and the phase opposition sum. For the
%latter, two statistical procedures are employed, a standard permutation
%test and a hybrid test using the mean and standard deviation of the null
%distribution obtained from the permutations to derive a z-score and,
%finally, a more precise p-value.
%The function operates on N-dimensional matrices of complex numbers.
%If data1/2 contain non-unit norm complex numbers, they will be normalized
%If data1/2 contain only real numbers, we assume that they represent phase
%angle (in radians).
%The function operates on the last (Nth) dimension, representing trials.
%Trials in data1 are associated with one task outcome, and in data2 with
%another task outcome. The first N-1 dimensions may represent
%anything, most likely time, frequency (obtained from any time-frequency
%transform) and electrodes/channels.
%
%For datasets with multiple subjects (and likely different trial numbers,
%run this function for each subject and combine the resulting maps of
%p-values using the function combine_pvalues (included in this folder).
%
%Inputs:
%   -data1: N-dimensional matrix of complex phase values or real angles (in
%   radians) for trial group 1
%   -data2: N-dimensional matrix of complex phase values or real angles (in
%   radians) for trial group 2
%   -[nperm]: number of permutations for non-parametric (POS) test
%   [default 1000]
%
%Outputs:
%   -p_POS: N-1-dimensional matrix of p-values obtained from the POS
%   measure and a comparison with a surrogate distribution of permutations
%   -p_zPOS: N-1-dimensional matrix of p-values obtained from the POS
%   measure and a zscore against the surrogate distribution of permutations

if nargin < 3
  nperm = 2;%1000;
end
N =1;

alldata = cat(N,data1,data2);

if max(abs(abs(alldata(:))-1))>0.00001 %non-unit norm
  fprintf('Complex numbers with non-unit norm, normalizing\n');
  alldata = alldata ./ abs(alldata);
  data1 = data1 ./ abs(data1);
  data2 = data2 ./ abs(data2);
end

%compute POS values (only sum, no baseline correction to facilitate
%surrogate comparisons)
fprintf('Computing POS and surrogates.');
apcd_A = ppc(data1);
apcd_B = ppc(data2);
appoi = squeeze(apcd_A + apcd_B);

%compute POS surrogates. To save memory, the distribution is not saved, but
%its mean and standard deviation, as well as the p_value count, are
%calculated iteratively.

% biased POS
mappoi = zeros(size(appoi)); %to store running mean
sappoi = zeros(size(appoi)); %to store running std
nappoi = zeros(size(appoi)); %to store running count

surrogate = zeros(1,nperm);
for s=1:nperm
  if ~rem(s-1,ceil(nperm/10)), fprintf('.'), end
  order = randperm(size(alldata,N));
  perm1 = order(1:size(data1,N));
  perm2 = order(size(data1,N)+1:end);
  
  %need to construct the index strings to pass to eval
  % biased appoi
  indexstring1 = 'alldata('; for i=1:N-1, indexstring1=[indexstring1 ':,']; end; indexstring1 = [indexstring1 'perm1)'];
  indexstring2 = 'alldata('; for i=1:N-1, indexstring2=[indexstring2 ':,']; end; indexstring2 = [indexstring2 'perm2)'];
  surrogateappoi = squeeze(ppc(eval(indexstring1)) + ppc(eval(indexstring2)));
  surrogate(1,s) = surrogateappoi;
  mappoi = mappoi+surrogateappoi;
  sappoi = sappoi+surrogateappoi.^2;
  nappoi = nappoi + double(surrogateappoi>appoi);
end

% biased appoi
mappoi = mappoi / nperm;
sappoi = sqrt((sappoi/nperm)-mappoi.^2);
nappoi(nappoi==0) = 0.5;

%compute probabilities from permutations
p = nappoi/nperm;
p_z = min(1-10^-20,max(10^-20,1-normcdf((appoi-mappoi)./sappoi)));
fprintf('\nbiased appoi done.\n')

function [c,a] = ppc(input)
n = size(input,1);
input  = input./abs(input);
outsum = nansum(input);
c = (outsum.*conj(outsum) - n)./(n*(n-1));
a = angle(mean(input));% angle at which the vector strength points at
% d2all = circ_dist2(angle(input),angle(input));
% mask = nonzeros(tril(d2all,-1));
% dx = abs(mean(exp(sqrt(-1).*mask)));
% ax = angle(mean(input));
% out = [];
% for j = 1:(n-1)
%   for k = (j+1):n
%     o1 = input(j,:).*conj(input(k,:));
% %   circ_r(nonzeros(tril(circ_dist2(angle(input),angle(input)),-1)))
%     out = [out;o1];
%   end
% end
% dx = abs(nansum(out,1)*2/(n*(n-1)));

% N = size(X,2);
% for j = 1:N-1
%     f = [];
%     for k = j+1:N
%         f(k) = cos(X(1,j))*cos(X(2,k)) + sin(X(1,j))*sin(X(2,k));
%     end
%     g(j) = sum(f);
% end
% h = (2*sum(g)) / (N*(N-1));