function [p, p_z, ppc_AB, ppc_A, ppc_B, surrogate] = pacoi(data1, data2, nperm)
% PACOI computes the Phase-Amplitude Coupling Opposition index (PACOi).
% Using the pairwise phase consistency (Vinck et al 2010 Neuroimage), it
% measures whether the PAC phase preference differs between two
% experimental conditions. The method is inspired from the Phase Opposition 
% Sum (VanRullen 2016 Front Neurosci).
%
% INPUT
%   data1       - N-dimensional matrix of complex phase values for trial
%               group 1. The trailing dimension must correspond to trials
%   data2       - N-dimensional matrix of complex phase values for trial
%               group 2. The trailing dimension must correspond to trials
%   nperm       - number of permutations to compute PACOi [default 1000]
%
% OUTPUT
%   p           - N-1-dimensional matrix of p-values obtained from the
%               PACOi measure and a comparison with a surrogate 
%               distribution of permutations
%   p_z         - N-1-dimensional matrix of p-values obtained from the
%               PACOi measure and a zscore against the surrogate
%               distribution of permutations
% 
% 08-Nov-2019 09:29:50 Diego Lozano-Soldevilla
%

if nargin < 3
  nperm = 1000;
end

N = ndims(data1);
if N~=ndims(data2)
  error('data1 and data2 should have the same number of dimensions\n');
end

sz1 = size(data1);
sz2 = size(data2);
if any(sz1(1:N-1)~=sz2(1:N-1))
  error('data1 and data2 should have the same 1:N-1 dimension size')
end

if any([isreal(data1),isreal(data2)])
  error('data1 and data2 shoulbe be complex\n');
end

alldata = cat(N,data1,data2);

if max(abs(abs(alldata(:))-1))>0.00001 %non-unit norm
  fprintf('Complex numbers with non-unit norm, normalizing\n');
  alldata = alldata ./ abs(alldata);
  data1 = data1 ./ abs(data1);
  data2 = data2 ./ abs(data2);
end

fprintf('Computing PACOi and surrogates\n');
ppc_A  = ppc(data1,N);
ppc_B  = ppc(data2,N);
ppc_AB = squeeze(ppc_A + ppc_B);

mppc_AB = zeros(size(ppc_AB)); %to store running mean
sppc_AB = zeros(size(ppc_AB)); %to store running std
nppc_AB = zeros(size(ppc_AB)); %to store running count

surrogate = nan([sz1(1:N-1) nperm]);
coldims = repmat(':,',1,N-1);
for s=1:nperm
  if ~rem(s-1,ceil(nperm/10)), fprintf('.'), end
  order = randperm(size(alldata,N));
  perm1 = order(1:size(data1,N));
  perm2 = order(size(data1,N)+1:end);
  indexstring1 = ['alldata(' coldims 'perm1)'];
  indexstring2 = ['alldata(' coldims 'perm2)'];
  surrogateppc_AB = ppc(eval(indexstring1),N) + ppc(eval(indexstring2),N);
  eval(['surrogate(' coldims 's) = surrogateppc_AB;']);
  mppc_AB = mppc_AB + surrogateppc_AB;
  sppc_AB = sppc_AB + surrogateppc_AB.^2;
  nppc_AB = nppc_AB + double(surrogateppc_AB > ppc_AB);
end

% biased ppc_AB
mppc_AB = mppc_AB / nperm;
sppc_AB = sqrt((sppc_AB/nperm)-mppc_AB.^2);
nppc_AB(nppc_AB==0) = 0.5;

%compute probabilities from permutations
p = nppc_AB/nperm;
p_z = min(1-10^-20,max(10^-20,1-normcdf((ppc_AB-mppc_AB)./sppc_AB)));
fprintf('\nbiased ppc_AB done.\n')

function c = ppc(input,dim)
n = size(input,dim);
input  = input./abs(input);
outsum = nansum(input,dim);
c = (outsum.*conj(outsum) - n)./(n*(n-1));

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
