%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Get MT projection of data...                                           %
%   % **** INPUT ****                                                        %
%   % data    - array of either lfp signals or of spike times                %          
%   % smp     - either the sampling rate if lfp or the freqs if spks         %  
%   % Ti      - interval for evaluation (start and end)                      %          
%   % W       - half bandwidth in Hz                                         %          
%   % pad     - fractional padding (eg 2 doubles grid only required for lfp) %
%   % kindx   - number of tapers used (use -1 for max allowed)               %
%   % flag    - 1 if lfp, 0 if spikes                                        %
%   % work     - precalculated tapers                                        %
%   %                                                                        %
%   % **** OUTPUT ****                                                       %
%   %                                                                        %
%   % J(trial,taper,f) - output matrix of projections                        %
%   % f                - frequencies                                         %
%   % hf               - high frequency limit (spikes)                       %
%   % work             - tapers for recall                                   %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[J,f,hf,work] = mtfft(data,smp,Ti,W,pad,kindx,flag,work);

% if in low spike, long time regime use explict representation of spikes...
% this can be a substantial memory saving and may also be required for 
% probing very high frequency structure (>1000Hz).  However because of the fft 
% it is quicker to put the spikes onto an fft grid under most circumstances
% note to enforce a preferred regime set regime here.  0 = use fft for spikes 
% 1 = use explict representation for spikes

regime = 0;
if flag == 0 & Ti(2)-Ti(1) > 60
  if length(find(data(1,:)>Ti(1) & data(1,:)<Ti(2)))/(Ti(2)-Ti(1)) < 5
    regime = 1; 
  end
end

% if it is quicker to treat spikes as a continuous process so put them on an fft grid...

if regime == 0 & flag == 0
  fq = smp;
  NT = length(data(:,1));
  datatmp = zeros(size(data));
  Np = zeros(NT,1);
  for n=1:NT
    indx = find(Ti(1)<data(n,:) & Ti(2)>data(n,:) & data(n,:) ~= 0);
    Np(n) = length(indx);
    datatmp(n,1:Np(n)) = data(n,indx);  
  end
  data = datatmp(:,1:max(Np));  
  smp = linspace(Ti(1),Ti(2)+0.001,fix(2000*(Ti(2)-Ti(1))));
  cdata = zeros(NT,length(smp));
  Dt = smp(2)-smp(1);  
  for n=1:NT
    indx = fix((data(n,1:Np(n))-Ti(1))/Dt) + 1;
    cdata(n,indx) = 1/Dt;
  end
  data = cdata;
end    
    
% check that work array paramters match previous ones...

if nargin == 8 & ~isempty(work)
 if ~((Ti(2)-Ti(1))-work.T) & ~(W-work.W) & ~(kindx-work.kindx) & ~(pad - work.pad) ...
                            & ~(flag-work.flag) & ~(regime-work.regime)
   E = work.E;
   H = work.H;
 else
   work = [];
 end
else
  work = [];
end

if isempty(work)
  work = struct('T',[],'W',[],'kindx',[],'pad',[],'flag',[],'regime',[],'E',[],'H',[],'t',[]);
  work.T = Ti(2)-Ti(1);
  work.W = W;
  work.kindx = kindx;
  work.pad = pad;  
  work.flag = flag;
  work.regime = regime;
end

verbose = 0;
hf = 0;
if flag == 1 | regime == 0 % continuous
  t = smp;
  smp_rate = 1/(t(2)-t(1));
  %N = fix((Ti(2)-Ti(1))*smp_rate);
  N = fix(((Ti(2)-Ti(1))*smp_rate)+.00000001);  % fix dumb matlab rounding error
  T = Ti(2)-Ti(1);
  tindx = find(t>Ti(1));
  if Ti(1) < t(1) | Ti(2) > t(length(t))
    disp('Range exceeds data')
    J = 0;
    return
  end

  tindx = tindx(1:N);
  
  DT = t(2)-t(1);
  NW = floor(N*W*DT);
  
  data = data(:,tindx);
  if kindx == -1 
    kindx = 2*NW -1;
  end
  if kindx < 1;error('Need a larger bandwidth or time interval');end
  
  if isempty(work.E)
    if verbose; disp('Setup tapers');end
    [E V] = dpss(N,NW,kindx);
    H = 0;
    work.E = E;
    work.t = linspace(0,Ti(2)-Ti(1),N);
  end
  
  if pad == 0; disp('Padding set to zero : reset to 2');pad = 2;end

  if pad > 0
    pad = fix(N*pad);
  else
    pad = -pad;
  end
  
  sz = size(data);
  dat = zeros(sz);
  if mod(pad,2)== 0
    jpad = fix(pad/2);
  else
    jpad = floor(pad/2)+1;
  end
  
  
  % new freq method, next pow 2
  %nf = 2^nextpow2(round(smp_rate));
  %f = [0:nf/2-1]*round(smp_rate)/nf;
  %pad = nf;
  %jpad = length(f);
  
  %keyboard
  
  % allocate space for data
  J = ones(sz(1),length(E(1,:)),jpad)*exp(i*eps);
  
  
  for k=1:length(E(1,:))
    for n=1:sz(1)
%      data(n,:) = data(n,:) - mean(data(n,:));
      data(n,:) = detrend(data(n,:));
      dat(n,:) = E(:,k)'.*data(n,:);
    end
    j = fft(dat,pad,2);
    J(:,k,:) = sqrt(2*DT)*j(:,1:jpad);
  end
  
  f = linspace(0,(pad-1)/(pad*DT),pad);
  f = f(1:jpad);
  
  %keyboard
  
  % if spikes extract the relevant freqs...
  if flag == 0 & regime == 0
    hf = Np/(Ti(2)-Ti(1));
    indx = zeros(1,length(fq));
    for n=1:length(fq)
      [mn indx(n)] = min((fq(n)-f).^2);  
    end
    indx = unique(indx);
    f = f(indx);
    J = J(:,:,indx)/sqrt(2);
  end
  
elseif flag == 0 & regime == 1                                % spikes

  f = smp;
  T = Ti(2)-Ti(1);
  Nf = length(f);
  sz = size(data);
  NT = sz(1);
  w = 2*pi*f; 
  datatmp = zeros(size(data));
  Np = zeros(1,NT);
  for n=1:NT
    indx = find(Ti(1)<data(n,:) & Ti(2)>data(n,:) & data(n,:) ~= 0);
    Np(n) = length(indx);
    datatmp(n,1:Np(n)) = data(n,indx) - Ti(1);  
  end
  data = datatmp(:,1:max(Np));  

  N = 10000;
  LL = (N-1)/T;
  NW = fix(2*W*T)/2;
  if kindx == -1
    kindx = floor(2*NW-1);
    if verbose; disp(['Used ' num2str(kindx) ' tapers']); end
  end
  if isempty(work.E)
    if verbose; disp('Setup tapers'); end
    [E V] = dpss(N,NW,kindx);
    E = E*sqrt(LL);
  end
  if isempty(work.H)
    if verbose; disp(['Transform tapers']); end
    t = linspace(0,T,N);
    H = ones(kindx,Nf)*exp(i*eps);
    for ff=1:Nf
      if f(ff) > 100*W; break; end
      ex = exp(-i*w(ff)*t);
      for k=1:kindx
        H(k,ff) = sum(E(:,k)'.*ex);
      end
    end
    H = H/LL;
  end
  work.E = E;
  work.H = H;
  
% enter main loop

  J = zeros(NT,kindx,Nf);
  hf = zeros(NT,1);
  for n=1:NT
    if verbose; disp(['Evaluate trial ' num2str(n)]); end
    for k=1:kindx
      if Np(n) > 0
        S = data(n,1:Np(n));
        MM = 1/LL;
        j_l = min(floor((S*LL) + 1),N-1);
        j_u = j_l + 1;
        t_l = (j_l-1)*MM;
        t_u = (j_u-1)*MM;  
	    h = (E(j_l,k)'.*(t_u-S)+ E(j_u,k)'.*(S-t_l))./(t_u-t_l);
        L = Np(n)/T;
        hf(n) = L;
        for ww = 1:Nf
          J(n,k,ww) = sum(h.*exp(-i*w(ww)*S)) - L*H(k,ww);
        end
      end
    end
  end
end







