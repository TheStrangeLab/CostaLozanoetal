%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Routine to remove line noise in spectra by fitting  %
%   % and removing a periodic signal.  The line if only   %
%   % removed if it is significant (Ftest)                %
%   % see Percival and Walden. Spectral analysis for      %
%   % physical applications fro details (this version for %
%   % continuous process data only                        %
%   %                                                     %
%   % INPUT                                               %
%   %                                                     %
%   % fq       - the frequencies at which the f-test is   %
%   %            performed (typically 60Hz + harmonics)   %
%   % data     - lfp data in form data(trials,sample)     %
%   % smp      - lfp sample times                         %
%   % plt      - 'y'|'n' (default 'y')                    % 
%   % T        - sub interval to analyze (default all)    %
%   % W        - bandwidth (Hz) (default 5)               %
%   % pad      - level of fft padding (default 2X)        %
%   % kindx    - tapers included -1 = all (default -1)    %
%   % err      - error bar type                           %
%   %            0 none                                   %
%   %            1 Asymptotic (default if trials < 5)     %
%   %            2 jackknife  (default if trials >= 5     %
%   % pvalue   - generate a 1-p value conf interval {0.05}%
%   % work     - precalculated taper functions            %
%   %                                                     %
%   % OUPUT                                               %
%   %                                                     %
%   % S        - residual spectrum                        %
%   % f        - frequencies                              %
%   % Up       - upper band of 95% conf interval          %
%   % Uq       - lower band of 95% conf interval          %
%   % p        - p value of fitted lines                  %
%   % work     - calculated tapers (use for speed)        %
%   % J        - set of tapered tranforms with lines      %
%   %            removed (useful for coherency etc)       %
%   %                                                     %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[S,f,Up,Uq,p,work,J] = ...
        cftest(fq,data,smp,plt,T,W,pad,kindx,err,pvalue,work)

verbose = 0;    
    
if isempty(fq); error('No frequencies requested');end    
if isempty(data); error('No data');end
if isempty(smp); error('I need the sample times');end
if nargin < 4; plt = 'r';end
if nargin < 5; T = [min(smp) max(smp)];end
if nargin < 6; W = 5; end
if nargin < 7; pad = 5; end
if nargin < 8; kindx = -1; end
if nargin < 9;if length(data(:,1)) > 4;err =2;else;err = 1;end;end
if nargin < 10; pvalue = 0.05;end
if nargin < 11; work = [];end

if isempty(plt);plt = 'r';end
if isempty(T);T = [min(smp) max(smp)];end
if isempty(W);W = 5;end
if isempty(pad);pad = 5;end
if isempty(kindx);kindx = -1;end
if isempty(err);if length(data(:,1)) > 4;err =2;else;err = 1;end;end
if isempty(pvalue);pvalue = 0.05;end

%  do the transform...

[J,f,hf,work] = mtfft(data,smp,T,W,pad,kindx,1,work);
J = J/sqrt(2);

% find the frequency closest to that tested...

for ff=1:length(fq)
  R = (fq(ff)-f).^2;
  findx(ff) = find(min(R)==R);
  fq(ff) = f(findx(ff));
end

% moved outside
E = work.E;
N = length(E(:,1));
if pad > 0
  pad = fix(N*pad);
else
  pad = -pad;
end
DT = smp(2)-smp(1);
Nf = length(f);
K = length(E(1,:));

% figure out the transforms of the tapers...
if isempty(work.H)
  H = ones(K,Nf)*exp(i*eps);
  for k=1:K
    jh = fft(E(:,k),pad,1)';
    H(k,1:Nf) = DT*conj(jh(1:Nf));
  end
  work.H = H;
else
  H = work.H;
end

%  now do the f-test...

NT = length(J(:,1,1));
for n=1:NT
  for nn=1:length(fq)
    if K < 3
      p(nn) = 100;
      if verbose
	disp(['Trial ' num2str(nn) ' has too few tapers '...
	      'for f-test'])
      end     
    else  
      for k=1:K
        H0(k) = H(k,1);
        J0(k) = J(n,k,findx(nn));
      end
      alpha = sum(H0.*conj(H0));
      c1(nn) = sqrt(DT)*sum(J0.*H0)/alpha;
      J_hat = c1(nn)*H0/sqrt(DT);
      F_num = abs(c1(nn)*conj(c1(nn))*alpha*(K-1));
      F_den(nn) = DT*abs(sum((J0-J_hat).*conj(J0-J_hat)));   
      F(nn) = F_num/F_den(nn);
      p(nn) = fix(10*(100*(1-fcdf(F(nn),2, ...
              2*(K-1)))))/10;
      Fd = fix(F(nn)*100)/100;
      fd = fix(fq(nn)*100)/100;
      if verbose
        disp(['trial = ' num2str(n) ...
            ' f = ' num2str(fd) ' Hz :     ' ...
            ' F = ' num2str(Fd)     ...
            ' p = ' num2str(p(nn)) '%'])
      end
    end
  end

% remove line from spectrum if significant at 20% level...
% this generous significance criteria is because it doesnt
% matter too much is a poor fit is removed but it is does
% if a true line is missed

  P = 20;                       % significance level
  Hf = ones(K,Nf)*exp(i*eps);
  if ~isempty(find(p<P))
    wf = 2*pi*fq;
    for nn=1:length(wf)
      if p(nn) < P  
	Hf(:,findx(nn):Nf) = H(:,1:(Nf-findx(nn)+1));
        Hf(:,1:(findx(nn)-1)) = ...
	    fliplr(conj(H(:,2:(findx(nn)))));
        J(n,:,:) = squeeze(J(n,:,:)) - c1(nn)*Hf/sqrt(DT);
      end
    end
  end
end
J = sqrt(2)*J;

%  average the estimates together (sum over tapers and trials)

sz = size(J);
S = zeros(sz(3),1);
for n=1:sz(1)
  for m=1:sz(2)
    S = S + abs(squeeze(J(n,m,:))).^2;
  end
end
dof = 2*sz(1)*sz(2);
S = S/(sz(1)*sz(2));
if err > 0
  if err == 1                     % chi2 CI
    pc = pvalue/2;
    q = 1-pc;
    Qp = chi2inv(pc,dof);
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
    Up = S.*exp(tinv((1-pvalue/2),N-1)*sig);
    Uq = S.*exp(-tinv((1-pvalue/2),N-1)*sig);   
  end
else
  Up = 0;
  Uq = 0;
end

% plot out the results...

if plt ~= 'n' 
  plot(f,S,[plt '-'])
  hold on
  if err > 0
    plot(f,Up,['g' '-'])
    plot(f,Uq,['g' '-'])
    hold off
  end
end









