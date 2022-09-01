%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Spectrum : calculates the MT spectrum                     %
%   %                                                           %
%   %                 ******** INPUT *********                  %
%   % Required...                                               %
%   %                                                           %
%   % data    - array of either lfp signals or of spike times   %          
%   %           should be of form data(trials,data)             %
%   % smp     - lfp times of samples                            %
%   %           spk frequencies of evaluation                   %
%   %                                                           %
%   % Optional...                                               %
%   %                                                           %
%   % Parameter         Example                 Default         %
%   %                                                           %
%   % plt (color or no plot) 'r'                 'n'            %
%   % interval T         [0.1 0.5]                all           %          
%   % W (bandwidth/2 Hz)   5                      5             %
%   % pad   factor of fft padding for lfp         2             %
%   % kindx number of tapers to be included       -1 (all)      %
%   % err   = 0 no errorbar                                     %
%   %       = 1 (asymptotic conf interval)                      %
%   %       = 2 (jackknife)                                     %
%   %       = 3 (approx finite size corrected theoretical)      %
%   % pvalue = p value of confidence interval                   %
%   % work data tapers (quicker if given)                       %
%   %                                                           %
%   % **** OUTPUT ****                                          %
%   % Spectrum S at frequencies (Hz) f                          %
%   % hf = high frequency limit (for point process spectra)     %
%   % Up,Uq is conf interval                                    %
%   % work are the data tapers (saves recalculation)            %   
%   % J  matrix of tapered transforms of the data               %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[S,f,hf,Up,Uq,work,J] = spectrum(data,smp,plt,T,W,pad,kindx,err,pvalue,work)

% figure out if it is lfp or spikes (flag = 0 indicates spikes) ...

if nargin < 2; error('data and sample times (lfp) or freqs of evaluation (spk) required');end

flag = 0;if length(find(diff(data(1,:))<0)) > 3;flag = 1;end

% set up defaults...

if nargin < 3;plt = 'n';end
if nargin < 4;if flag;T= [min(smp) max(smp)];else;T= [min(data(:,1)) max(max(data))];end;end
if nargin < 5;  W = 5; end
if nargin < 6; pad = 2; end
if nargin < 7; kindx = -1; end
if nargin < 8;if length(data(:,1)) > 4;err = 2;else;err = 1;end;end
if nargin < 9; pvalue = 0.05;end
if nargin < 10;work = [];end

if isempty(plt);plt = 'n';end
if isempty(T);if flag;T= [min(smp) max(smp)];else;T= [min(data(:,1)) max(max(data))];end;end
if isempty(W);  W = 5; end
if isempty(pad); pad = 2; end
if isempty(kindx); kindx = -1; end
if isempty(err);if length(data(:,1)) > 4;err = 2;else;err = 1;end;end
if isempty(pvalue); pvalue = 0.05;end

%  do the transform...

[J,f,hf,work] = mtfft(data,smp,T,W,pad,kindx,flag,work);
S = squeeze(mean(mean(abs(J).^2,1),2));

%  finite size correction to errorbars...

sz = size(J);
dof = 2*sz(1)*sz(2);

if dof == 0; error(['No degrees of freedom! Increase bandwidth or' ...
	' time interval']);end

if flag == 0 & err == 3
  totspk = 0;
  for n=1:sz(1)
    indx = find(T(1)<data(n,:) & T(2)>data(n,:) & data(n,:) ~= 0);
    totspk = totspk + length(indx);
  end
  if totspk > 0; dof = fix(1/(1/dof + 1/(2*totspk))); end 
end

if err > 0
  if err == 1 | err == 3                    % chi2 CI
    p = pvalue/2;
    q = 1-p;
    Qp = chi2inv(p,dof);
    Qq = chi2inv(q,dof);
    Up = dof.*S./Qp;
    Uq = dof.*S./Qq;
  elseif err == 2
    N = sz(1)*sz(2);
    mnlS = 0;
    sqlS = 0;
    for n=1:sz(1)
      for m=1:sz(2)
        lSmn = log((N*S - abs(squeeze(J(n,m,:))).^2)/(N-1));
        mnlS = mnlS + lSmn;
	    sqlS = sqlS + lSmn.^2;
      end
    end
    mnlS = mnlS/N;
    sqlS = sqlS/N;
    sig = sqrt((N-1)*(sqlS - mnlS.^2));
    Up = S.*exp(tinv((1-0.025),N-1)*sig);
    Uq = S.*exp(-tinv((1-0.025),N-1)*sig);    
  end
else
  Up = 0;
  Uq = 0;
end

if plt ~= 'n' 
  plot(f,S,[plt '-'])
  hold on
  if err > 0
    plot(f,Up,['g' '-'])
    plot(f,Uq,['g' '-'])
    hold off
  end
  if flag == 0
    line(get(gca,'xlim'),mean(hf)*[1 1],'color','b')
    set(gca,'ylim',[0 1.25*max(S)])
    set(gca,'xlim',[min(f) max(f)])
  end
end
























