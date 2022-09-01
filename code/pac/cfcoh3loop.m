function [MI,z,k2sel] = cfcoh3loop(raw1,raw2,raw3,flow,flow_bw,fhigh,fhigh_bw,position,method,catrpt)
%
% CFCOH3LOOP compute the pair-wise phase consistency (Vinck et al 2010
% Neuroimage) between the phase difference taken between two time courses
% (channel 1 - channel 2) as a function of a third phase time course (channel 3).
%
% INPUT
%  	raw1        - unfiltered signal (channed 1) in raw format (see ft_datatype_raw)
%               that will be filtered by fhigh frequencies
%  	raw2        - unfiltered signal (channel 2) in raw format (see ft_datatype_raw)
%               that will be filtered by fhigh frequencies
%   raw3        - unfiltered signal (channel 3) in raw format (see ft_datatype_raw)
%               whose phase will be binned by flow frequencies
%   flow        - low frequency axis (1xN vector)
%   flow_bw     - low frequency bandwidth such that [flow(x)-flow_bw/2 flow(x)+flow_bw/2];  (number)
%   flow        - high frequency axis (1xN vector)
%   fhigh_bw    - high frequency bandwidth such that [fhigh(y)-fhigh_bw/2 fhigh(y)+fhigh_bw/2]; (number)
%   position    - low frequency phase binning (number)
%   method      - estimate the phase difference between the filtered raw1
%               and raw2 time courses ('phase') or the phase difference
%               between the power envelopes of the filtered raw1 and raw2
%               ('powenv')
%   catrpt      - 0 or 1. Concatenate single trials after filtering to
%               increase the number of samples to compute PPC
%
% OUTPUT
%   MI          - modulation index computed over the raw1-raw2 as a function of raw3
%   z           - complex representation of pairwise phase consistency
%               between raw1-raw2 as a function of raw3
%   k2sel       - trials included in the analysis
%
%  20-Dec-2020 01:05:33, Diego Lozano-Soldevilla
%

dat1 = raw2data(raw1,'rpt_chan_time');
dat1 = squeeze(dat1.trial);

fsample = raw1.fsample;

dat2 = raw2data(raw2,'rpt_chan_time');
dat2 = squeeze(dat2.trial);

dat3 = raw2data(raw3,'rpt_chan_time');
dat3 = squeeze(dat3.trial);

clear raw1 raw2 raw3;

if size(dat1)~=size(dat2) | size(dat1)~=size(dat3) | size(dat2)~=size(dat3)
  error('dat with different dimensions')
end
[nrpt,nsmp] = size(dat1);

nbin=length(position);
winsize = 2*pi/nbin;
[phs_mcat, phs1cat, phs2cat]=deal([]);% in case

if catrpt==0
  ampbin = nan(nrpt,size(flow,2),size(fhigh,2),nbin);
else
  ampbin = nan(size(flow,2),size(fhigh,2),nbin);
end
for x = 1:size(flow,2)
  Fbp_x = [flow(x)-flow_bw/2 flow(x)+flow_bw/2];
  N_x = 3*fix(fsample/Fbp_x(1));
  
  for y = 1:size(fhigh,2)
    Fbp_y = [fhigh(y)-fhigh_bw/2 fhigh(y)+fhigh_bw/2];
    N_y = 3*fix(fsample/Fbp_y(1));
    
    for k = 1:nrpt
      filt_x = ft_preproc_bandpassfilter(zscore(dat3(k,:)),fsample,Fbp_x,N_x,'fir','twopass');
      phs_m = ft_preproc_hilbert(filt_x,'angle');
      filt_y1 = ft_preproc_bandpassfilter(zscore(dat1(k,:)),fsample,Fbp_y,N_y,'fir','twopass');
      filt_y2 = ft_preproc_bandpassfilter(zscore(dat2(k,:)),fsample,Fbp_y,N_y,'fir','twopass');
      
      switch method
        case 'phase'
          phs1 = ft_preproc_hilbert(filt_y1,'angle');
          phs2 = ft_preproc_hilbert(filt_y2,'angle');
        case 'powenv'
          pow1 = ft_preproc_hilbert(filt_y1,'abs').^2;
          pow2 = ft_preproc_hilbert(filt_y2,'abs').^2;
          phs1 = ft_preproc_hilbert(pow1,'angle');
          phs2 = ft_preproc_hilbert(pow2,'angle');
      end
      
      if catrpt == 0
        tmp = nan(1,nbin);
        for j=1:nbin
          phs_mask = find(phs_m <  position(j)+winsize & phs_m >=  position(j));
          tmp(1,j) = ppc(exp(sqrt(-1)*(phs1(phs_mask) - phs2(phs_mask))),2);
        end
        phs_mask=[];
        ampbin(k,x,y,:)=tmp;
      else
        phs_mcat = cat(2,phs_mcat,phs_m);
        phs1cat = cat(2,phs1cat,phs1);
        phs2cat = cat(2,phs2cat,phs2);
      end
    end
    
    if catrpt==1
      tmp = nan(1,nbin);
      for j=1:nbin
        phs_mask = find(phs_mcat <  position(j)+winsize & phs_mcat >=  position(j));
        tmp(1,j) = ppc(exp(sqrt(-1)*(phs1cat(phs_mask) - phs2cat(phs_mask))),2);
      end
      phs_mask=[];
      ampbin(x,y,:)=tmp;
    end
    
  end
end

if catrpt==1
  z = [];k2sel=[];% as trials are concatenated I leave this empty
  MeanAmp = ampbin;
  nAmp = MeanAmp./sum(MeanAmp,3);
  nAmp(nAmp<0)=nan; % sometimes the ppc gives small but negative values that I ignore
  % https://doi.org/10.1016/j.neuroimage.2010.01.073
  %  Another potential disadvantage is that negative values are possible when
  %  using the normalized version of the PPC (whose expected value runs from
  %  0 to 1). This is a consequence of the unbiasedness of the PPC. In
  %  general, if a statistic is bounded from below by zero, then its expected
  %  value can only be zero then the variance is zero.
  Hp = -nansum(nAmp.*log(nAmp),3);
  Dkl = log(nbin)-Hp;
  MI = Dkl/log(nbin);
else
  k2sel = find(sum(isnan(ampbin(:,:)),2)==0);
  clear tmp;
  ampbin = ampbin(k2sel,:,:,:);
  siz = size(ampbin);
  z = nan(siz);
  for f1=1:siz(2)
    for f2=1:siz(3)
      for k=1:siz(1)
        z(k,f1,f2,:) = squeeze(ampbin(k,f1,f2,:))'.*exp(sqrt(-1).*position);
      end
    end
  end
  MeanAmp = squeeze(mean(ampbin,1));
  nAmp = MeanAmp./sum(MeanAmp,3);
  nAmp(nAmp<0)=nan;
  Hp = -nansum(nAmp.*log(nAmp),3);
  Dkl = log(nbin)-Hp;
  MI = Dkl/log(nbin);
end

function c = ppc(input,dim)
n = size(input,dim);
input  = input./abs(input);
outsum = nansum(input,dim);
c = (outsum.*conj(outsum) - n)./(n*(n-1));

function [data] = raw2data(data, dimord)

% RAW2DATA is a helper function that converts raw data to various types of
% averages. This function is used to apply the analysis steps that were
% written for use on preprocessed data also on averaged data.
%
% This function is the counterpart of DATA2RAW and is used in MEGREALIGN, MEGPLANAR, MEGREPAIR

% Copyright (C) 2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if isempty(dimord)
  % no conversion is needed
  return;
end

switch dimord
  case 'chan_time'
    fprintf('converting single trial back to average\n');
    data.avg    = data.trial{1};
    data.time   = data.time{1};
    data        = rmfield(data, 'trial');
    data.dimord = dimord;
    
  case 'rpt_chan_time'
    fprintf('converting raw trials to timelocked trials\n');
    ntrial = length(data.trial);
    nchan  = size(data.trial{1},1);
    ntime  = size(data.trial{1},2);
    tmptrial = zeros(ntrial, nchan, ntime);
    for i=1:ntrial
      tmptrial(i,:,:) = data.trial{i};
    end
    data = rmfield(data, 'trial');
    data.trial  = tmptrial;
    data.avg    = reshape(mean(tmptrial, 1), nchan, ntime);         % recompute the average
    data.var    = reshape(std(tmptrial, [], 1).^2, nchan, ntime);   % recompute the variance
    data.time   = data.time{1};                                     % the time axes of all trials are the same
    data.dimord = dimord;
    
  case 'subj_chan_time'
    fprintf('converting raw trials to individual subject averages\n');
    nsubj  = length(data.trial);
    nchan  = size(data.trial{1},1);
    ntime  = size(data.trial{1},2);
    tmptrial = zeros(nsubj, nchan, ntime);
    for i=1:nsubj
      tmptrial(i,:,:) = data.trial{i};
    end
    data = rmfield(data, 'trial');
    data.individual = tmptrial;
    data.avg        = reshape(mean(tmptrial, 1), nchan, ntime);         % recompute the average
    data.var        = reshape(std(tmptrial, [], 1).^2, nchan, ntime);   % recompute the variance
    data.time       = data.time{1};                                     % the time axes of all trials are the same
    data.dimord     = dimord;
    
  otherwise
    ft_warning('unrecognized dimord');
end

