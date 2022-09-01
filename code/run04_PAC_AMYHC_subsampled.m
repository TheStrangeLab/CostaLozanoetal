% This script computes the phase to amplitude coupling (PAC) between amygdala
% phase and hippocampus amplitude controlling the number of trials. Supplementary Fig.11,13 

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20190203\'
ft_defaults

clear all,close all,clc
rng(144);

mlist={'s6 right', 's13', 's15', 's16', 's16 right', 's25', 's32', 's10', 's33'};%;

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','eeg_toolbox_v1_3_2'))% kahana eeg toolbox
addpath(fullfile(oripath,'Costalozanoetal','code','utils','MatlabImportExport_v6.0.0'));
addpath(fullfile(oripath,'Costalozanoetal','code','utils','circstat-matlab'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','teg_repeated_measures_ANOVA'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','violinPlot'))
addpath(fullfile(oripath,'Costalozanoetal','code','pac'))

load(fullfile(oripath,'data2github','DataCleanAH_RKF.mat'))
% column meaning: eR_hc eKF_hc nR_hc nKF_hc
load(fullfile(oripath,'data2github','trials_patients.mat'))

toi = [0.4102 1.1006];
dataeRclean = prepare_iapsKF_9patients(dataeR_AH);clear dataeR_AH;
dataeFclean = prepare_iapsKF_9patients(dataeF_AH);clear dataeF_AH;
dataeKclean = prepare_iapsKF_9patients(dataeK_AH);clear dataeK_AH;
datanRclean = prepare_iapsKF_9patients(datanR_AH);clear datanR_AH;
datanFclean = prepare_iapsKF_9patients(datanF_AH);clear datanF_AH;
datanKclean = prepare_iapsKF_9patients(datanK_AH);clear datanK_AH;

for subj = 4:8
  cfg = [];
  cfg.channel = 'A1';
  cfg.latency = toi;
  eR_amy{subj} = ft_selectdata(cfg, dataeRclean{subj});
  eF_amy{subj} = ft_selectdata(cfg, dataeFclean{subj});
  eK_amy{subj} = ft_selectdata(cfg, dataeKclean{subj});
  nR_amy{subj} = ft_selectdata(cfg, datanRclean{subj});
  nF_amy{subj} = ft_selectdata(cfg, datanFclean{subj});
  nK_amy{subj} = ft_selectdata(cfg, datanKclean{subj});
end

for subj = [1:3,9]
  cfg = [];
  cfg.channel = 'A2';
  cfg.latency = toi;
  eR_amy{subj} = ft_selectdata(cfg, dataeRclean{subj});
  eF_amy{subj} = ft_selectdata(cfg, dataeFclean{subj});
  eK_amy{subj} = ft_selectdata(cfg, dataeKclean{subj});
  nR_amy{subj} = ft_selectdata(cfg, datanRclean{subj});
  nF_amy{subj} = ft_selectdata(cfg, datanFclean{subj});
  nK_amy{subj} = ft_selectdata(cfg, datanKclean{subj});
end

for subj = [1,3,6:9]
  cfg = [];
  cfg.channel = 'HC1';
  cfg.latency = toi;
  eR_hc{subj} = ft_selectdata(cfg, dataeRclean{subj});
  eF_hc{subj} = ft_selectdata(cfg, dataeFclean{subj});
  eK_hc{subj} = ft_selectdata(cfg, dataeKclean{subj});
  nR_hc{subj} = ft_selectdata(cfg, datanRclean{subj});
  nF_hc{subj} = ft_selectdata(cfg, datanFclean{subj});
  nK_hc{subj} = ft_selectdata(cfg, datanKclean{subj});
end

for subj = [2,4,5]
  cfg = [];
  cfg.channel = 'HC2';
  cfg.latency = toi;
  eR_hc{subj} = ft_selectdata(cfg, dataeRclean{subj});
  eF_hc{subj} = ft_selectdata(cfg, dataeFclean{subj});
  eK_hc{subj} = ft_selectdata(cfg, dataeKclean{subj});
  nR_hc{subj} = ft_selectdata(cfg, datanRclean{subj});
  nF_hc{subj} = ft_selectdata(cfg, datanFclean{subj});
  nK_hc{subj} = ft_selectdata(cfg, datanKclean{subj});
end

cfg= [];
for subj = 1:9
  e_amy{subj} = ft_appenddata(cfg, eR_amy{subj}, eK_amy{subj}, eF_amy{subj});
  n_amy{subj} = ft_appenddata(cfg, nR_amy{subj}, nK_amy{subj}, nF_amy{subj});
  e_hc{subj}  = ft_appenddata(cfg, eR_hc{subj},  eK_hc{subj},  eF_hc{subj});
  n_hc{subj}  = ft_appenddata(cfg, nR_hc{subj},  nK_hc{subj},  nF_hc{subj});
end

for subj = 1:9
  cfg= [];
  eKF_amy{subj} = ft_appenddata(cfg, eK_amy{subj}, eF_amy{subj});
  nKF_amy{subj} = ft_appenddata(cfg, nK_amy{subj}, nF_amy{subj});
  eKF_hc{subj}  = ft_appenddata(cfg, eK_hc{subj}, eF_hc{subj});
  nKF_hc{subj}  = ft_appenddata(cfg, nK_hc{subj}, nF_hc{subj});
end
clear eK_amy eF_amy nK_amy nF_amy eK_hc eF_hc nK_hc nF_hc

cfg= [];
for subj = 1:9
  r_amy{subj} = ft_appenddata(cfg, eR_amy{subj},  nR_amy{subj});
  f_amy{subj} = ft_appenddata(cfg, eKF_amy{subj}, nKF_amy{subj});
  r_hc{subj}  = ft_appenddata(cfg, eR_hc{subj},   nR_hc{subj});
  f_hc{subj}  = ft_appenddata(cfg, eKF_hc{subj},  nKF_hc{subj});
end


srate = 500;
dt = 1/srate;

nbin = 20; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin
  position(j) = -pi+(j-1)*winsize;
end

flow  = 5:1:20;
fhigh = 40:5:120;
flow_bw  = 4;
fhigh_bw = 30;

% emotion effect variables
eMIx = nan(size(flow,2),size(fhigh,2),9);
eAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);

nMIx = nan(size(flow,2),size(fhigh,2),9);
nAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);

% memory effect variables
rMIx = nan(size(flow,2),size(fhigh,2),9);
rAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);

fMIx = nan(size(flow,2),size(fhigh,2),9);
fAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);

% interaction variables
eRMIx = nan(size(flow,2),size(fhigh,2),9);
eRAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);

eFMIx = nan(size(flow,2),size(fhigh,2),9);
eFAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);

nRMIx = nan(size(flow,2),size(fhigh,2),9);
nRAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);

nFMIx = nan(size(flow,2),size(fhigh,2),9);
nFAm  = nan(size(flow,2),size(fhigh,2),size(position,2),9);


nresamp = 20;
% trials: patients(9,1) x conditions (eR_hc, eKF_hc, nR_hc, nKF_hc)
trl=[];

for s=1:size(eR_hc,2)
  aeR = raw2data(eR_amy{1,s},'rpt_chan_time');
  aeR = squeeze(aeR.trial);
  heR = raw2data(eR_hc{1,s}, 'rpt_chan_time');
  heR = squeeze(heR.trial);
  [mi,meanamp1] = ModIndex_v2_ft_subsamp(aeR,heR,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel1 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,1)=size(k2sel1,1);
  MeanAmp = (mean(meanamp1(:,:,:,k2sel1),4));
  eRAm(:,:,:,s) = MeanAmp;
  eRMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];

  aeKF = raw2data(eKF_amy{1,s},'rpt_chan_time');
  aeKF = squeeze(aeKF.trial);
  heKF = raw2data(eKF_hc{1,s}, 'rpt_chan_time');
  heKF = squeeze(heKF.trial);

  [mi,meanamp2] = ModIndex_v2_ft_subsamp(aeKF,heKF,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel2 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,2)=size(k2sel2,1);
  MeanAmp = (mean(meanamp2(:,:,:,k2sel2),4));
  eFAm(:,:,:,s) = MeanAmp;
  eFMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];

  anR = raw2data(nR_amy{1,s},'rpt_chan_time');
  anR = squeeze(anR.trial);
  hnR = raw2data(nR_hc{1,s}, 'rpt_chan_time');
  hnR = squeeze(hnR.trial);

  [mi,meanamp3] = ModIndex_v2_ft_subsamp(anR,hnR,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel3 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,3)=size(k2sel3,1);
  MeanAmp = (mean(meanamp3(:,:,:,k2sel3),4));
  nRAm(:,:,:,s) = MeanAmp;
  nRMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];

  anKF = raw2data(nKF_amy{1,s},'rpt_chan_time');
  anKF = squeeze(anKF.trial);
  hnKF = raw2data(nKF_hc{1,s}, 'rpt_chan_time');
  hnKF = squeeze(hnKF.trial);

  [mi,meanamp4] = ModIndex_v2_ft_subsamp(anKF,hnKF,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel4 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,4)=size(k2sel4,1);
  MeanAmp = (mean(meanamp4(:,:,:,k2sel4),4));
  nFAm(:,:,:,s) = MeanAmp;
  nFMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];

  % main effect of emotion to include
  ae = raw2data(e_amy{1,s},'rpt_chan_time');
  ae = squeeze(ae.trial);
  he = raw2data(e_hc{1,s}, 'rpt_chan_time');
  he = squeeze(he.trial);
  [mi,meanamp5] = ModIndex_v2_ft_subsamp(ae,he,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel5 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,5)=size(k2sel5,1);
  MeanAmp = (mean(meanamp5(:,:,:,k2sel5),4));
  eAm(:,:,:,s) = MeanAmp;
  eMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];

  an = raw2data(n_amy{1,s},'rpt_chan_time');
  an = squeeze(an.trial);
  hn = raw2data(n_hc{1,s}, 'rpt_chan_time');
  hn = squeeze(hn.trial);
  [mi,meanamp6] = ModIndex_v2_ft_subsamp(an,hn,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel6 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,6)=size(k2sel6,1);
  MeanAmp = (mean(meanamp6(:,:,:,k2sel6),4));
  nAm(:,:,:,s) = MeanAmp;
  nMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];

  % memory effect variables
  ar = raw2data(r_amy{1,s},'rpt_chan_time');
  ar = squeeze(ar.trial);
  hr = raw2data(r_hc{1,s}, 'rpt_chan_time');
  hr = squeeze(hr.trial);
  [mi,meanamp7] = ModIndex_v2_ft_subsamp(ar,hr,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel7 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,7)=size(k2sel7,1);
  MeanAmp = (mean(meanamp7(:,:,:,k2sel7),4));
  rAm(:,:,:,s) = MeanAmp;
  rMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];

  af = raw2data(f_amy{1,s},'rpt_chan_time');
  af = squeeze(af.trial);
  hf = raw2data(f_hc{1,s}, 'rpt_chan_time');
  hf = squeeze(hf.trial);
  [mi,meanamp8] = ModIndex_v2_ft_subsamp(af,hf,srate,flow,flow_bw,fhigh,fhigh_bw,position,min(trials(s,:)),nresamp);

  mi2 = permute(mi,[3 1 2]);
  k2sel8 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,8)=size(k2sel8,1);
  MeanAmp = (mean(meanamp8(:,:,:,k2sel8),4));
  fAm(:,:,:,s) = MeanAmp;
  fMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];MeanAmp=[];
  clear meanamp*
end


%%
eRfreq=[];
eRfreq.label={'chan1'};
eRfreq.freq = fhigh;
eRfreq.time = flow;
eRfreq.powspctrm(:,1,:,:) = permute(eRMIx,[3 2 1]);
eRfreq.dimord = 'subj_chan_freq_time';

eFfreq = eRfreq;
eFfreq.powspctrm(:,1,:,:) = permute(eFMIx,[3 2 1]);

nRfreq = eRfreq;
nRfreq.powspctrm(:,1,:,:) = permute(nRMIx,[3 2 1]);

nFfreq = eRfreq;
nFfreq.powspctrm(:,1,:,:) = permute(nFMIx,[3 2 1]);

% emotion main effect
efreq = eRfreq;
efreq.powspctrm(:,1,:,:) = permute(eMIx,[3 2 1]);

nfreq = eRfreq;
nfreq.powspctrm(:,1,:,:) = permute(nMIx,[3 2 1]);

% memory main effect
rfreq = eRfreq;
rfreq.powspctrm(:,1,:,:) = permute(rMIx,[3 2 1]);

ffreq = eRfreq;
ffreq.powspctrm(:,1,:,:) = permute(fMIx,[3 2 1]);

eRfreq.powspctrm(8,:,:,:)=[];
eFfreq.powspctrm(8,:,:,:)=[];
nRfreq.powspctrm(8,:,:,:)=[];
nFfreq.powspctrm(8,:,:,:)=[];

%%
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract';
emo_d = ft_math(cfg, eRfreq, eFfreq);
neu_d = ft_math(cfg, nRfreq, nFfreq);


hfoi=[40 120];
lfoi=[5 20];

cfg = [];
cfg.method           = 'ft_statistics_montecarlo'; %'analytic';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusterstatistic = 'maxsum';
cfg.latency          = lfoi;
cfg.avgovertime      = 'no';
cfg.frequency        = hfoi;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;
cfg.neighbours = [];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,8) ones(1,8).*2; [1:8] [1:8]];
stat = ft_freqstatistics(cfg,emo_d,neu_d);

cfg.design = [ones(1,9) ones(1,9).*2; [1:9] [1:9]];
statEmo9 = ft_freqstatistics(cfg,efreq,nfreq);

efreq8 = efreq;
nfreq8 = nfreq;
efreq8.powspctrm(8,:,:,:)=[];
nfreq8.powspctrm(8,:,:,:)=[];

cfg.design = [ones(1,8) ones(1,8).*2; [1:8] [1:8]];
statEmo8 = ft_freqstatistics(cfg,efreq8,nfreq8);

cfg.design = [ones(1,9) ones(1,9).*2; [1:9] [1:9]];
statMem = ft_freqstatistics(cfg,rfreq,ffreq);

rfreq8 = rfreq;
ffreq8 = ffreq;
rfreq8.powspctrm(8,:,:,:)=[];
ffreq8.powspctrm(8,:,:,:)=[];

cfg.design = [ones(1,8) ones(1,8).*2; [1:8] [1:8]];
statMem8 = ft_freqstatistics(cfg,rfreq8,ffreq8);

diffmem9 = rfreq;
diffmem9.powspctrm = rfreq.powspctrm-ffreq.powspctrm;

diffmem8 = rfreq8;
diffmem8.powspctrm = rfreq8.powspctrm-ffreq8.powspctrm;

diff = emo_d;
diff.powspctrm = emo_d.powspctrm-neu_d.powspctrm;

diff9 = efreq;
diff9.powspctrm = efreq.powspctrm-nfreq.powspctrm;

diff8 = efreq8;
diff8.powspctrm = efreq8.powspctrm-nfreq8.powspctrm;

xemo_d = ft_freqdescriptives([],emo_d);
xneu_d = ft_freqdescriptives([],neu_d);
xdiff  = ft_freqdescriptives([],diff);

xefreq = ft_freqdescriptives([],efreq);
xnfreq = ft_freqdescriptives([],nfreq);
xdiff9 = ft_freqdescriptives([],diff9);

xefreq8 = ft_freqdescriptives([],efreq8);
xnfreq8 = ft_freqdescriptives([],nfreq8);
xdiff8 = ft_freqdescriptives([],diff8);

xrfreq    = ft_freqdescriptives([],rfreq);
xffreq    = ft_freqdescriptives([],ffreq);
xdiffmem9 = ft_freqdescriptives([],diffmem9);

xrfreq8    = ft_freqdescriptives([],rfreq8);
xffreq8    = ft_freqdescriptives([],ffreq8);
xdiffmem8 = ft_freqdescriptives([],diffmem8);

cfg=[];
cfg.latency = lfoi;
cfg.frequency = hfoi;
xemo_d = ft_selectdata(cfg,xemo_d);
xemo_d.mask = (stat.negclusterslabelmat==1);
xneu_d = ft_selectdata(cfg,xneu_d);
xneu_d.mask = (stat.negclusterslabelmat==1);
xdiff = ft_selectdata(cfg,xdiff);
xdiff.mask = (stat.negclusterslabelmat==1);

xefreq = ft_selectdata(cfg,xefreq);
xefreq.mask = (statEmo9.posclusterslabelmat==1);
xnfreq = ft_selectdata(cfg,xnfreq);
xnfreq.mask = (statEmo9.posclusterslabelmat==1);
xdiff9 = ft_selectdata(cfg,xdiff9);
xdiff9.mask = (statEmo9.posclusterslabelmat==1);

xefreq8 = ft_selectdata(cfg,xefreq8);
xefreq8.mask = (statEmo8.posclusterslabelmat==1);
xnfreq8 = ft_selectdata(cfg,xnfreq8);
xnfreq8.mask = (statEmo8.posclusterslabelmat==1);
xdiff8 = ft_selectdata(cfg,xdiff8);
xdiff8.mask = (statEmo8.posclusterslabelmat==1);

% memory main effect
xrfreq = ft_selectdata(cfg,xrfreq);
xrfreq.mask = (statMem.posclusterslabelmat==1);
xffreq = ft_selectdata(cfg,xffreq);
xffreq.mask = (statMem.posclusterslabelmat==1);
xdiffmem9 = ft_selectdata(cfg,xdiffmem9);
xdiffmem9.mask = (statMem.posclusterslabelmat==1);

xrfreq8 = ft_selectdata(cfg,xrfreq8);
xrfreq8.mask = (statMem8.posclusterslabelmat==1);
xffreq8 = ft_selectdata(cfg,xffreq8);
xffreq8.mask = (statMem8.posclusterslabelmat==1);
xdiffmem8 = ft_selectdata(cfg,xdiffmem8);
xdiffmem8.mask = (statMem8.posclusterslabelmat==1);

%%
cfg = [];
cfg.figure = 'gcf';
cfg.parameter = 'powspctrm';
cfg.maskparameter='mask';
cfg.maskstyle = 'opacity';
cfg.maskalpha = 1;
cfg.colormap = 'hot';
cfg.zlim = [0 1e-3];

h1 = Figure(1,'size',[120   120],'fontSizeScreen',10); %set(h1, 'Position', [68     4   979   657]);
subplot(331);ft_singleplotTFR(cfg,xemo_d);
tit=title('emo R - emo F');set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
h = colorbar;
h.FontSize = 8;
ylabel(h, 'MI')
subplot(332);ft_singleplotTFR(cfg,xneu_d);
tit=title('neu R - neu F');set(findobj(tit,'type','text'),'FontSize',10);
h=colorbar;
h.FontSize = 8;
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
ylabel(h, 'MI')

cfg.zlim = [-5 5];
cfg.parameter = 'stat';
cfg.maskalpha = 0.3;
stat.mask = (stat.negclusterslabelmat==1);
subplot(333);ft_singleplotTFR(cfg,stat);
tit=title(['Interact; pval=' num2str(stat.negclusters(1).prob)]);set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
h = colorbar;
h.FontSize=6;
ylabel(h, 'T-score ')
colormap(gca,'default')

% emotion main effect with 8 and 9 patients
cfg = [];
cfg.figure = 'gcf';
cfg.parameter = 'powspctrm';
cfg.maskparameter='mask';
cfg.maskstyle = 'opacity';
cfg.maskalpha = 1;
cfg.zlim = [0 1e-3];

xefreq.powspctrm(~xefreq.mask)=0;
subplot(3,3,4);ft_singleplotTFR(cfg,xefreq);
tit=title('Emotional');set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
colormap(gca,'hot')
h = colorbar;
h.FontSize = 8;
ylabel(h, 'MI')

xnfreq.powspctrm(~xnfreq.mask)=0;
subplot(3,3,5);ft_singleplotTFR(cfg,xnfreq);
tit=title('Neutral');set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
colormap(gca,'hot')
h = colorbar;
h.FontSize = 8;
ylabel(h, 'MI')

cfg.parameter = 'stat';
cfg.zlim = [-5 5];
cfg.maskalpha = 0.3;
statEmo9.mask = (statEmo9.posclusterslabelmat==1);
subplot(3,3,6);ft_singleplotTFR(cfg,statEmo9);
tit=title(['pval=' num2str(statEmo9.posclusters(1).prob)]);set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)')
ylabel('HC Amplitude (Hz)','FontSize',10);
colormap(gca,'default')
h = colorbar;
h.FontSize = 8;
ylabel(h, 'T-score')


% memory main effect with 8 and 9 patients
cfg = [];
cfg.figure = 'gcf';
cfg.parameter = 'powspctrm';
cfg.maskparameter='mask';
cfg.maskstyle = 'opacity';
cfg.maskalpha = 1;
cfg.zlim = [0 3e-4];
xrfreq.powspctrm(~xrfreq.mask)=0;

subplot(3,3,7);ft_singleplotTFR(cfg,xrfreq);
tit=title('Remembered');set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
colormap(gca,'hot')
h = colorbar;
h.FontSize = 8;
ylabel(h, 'MI')
xffreq.powspctrm(~xffreq.mask)=0;

subplot(3,3,8);ft_singleplotTFR(cfg,xffreq);
tit=title('Forgotten');set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
colormap(gca,'hot')
h = colorbar;
h.FontSize = 8;
ylabel(h, 'MI')

cfg.parameter = 'stat';
cfg.maskalpha = 0.3;
cfg.zlim = [-5 5];
statMem.mask = (statMem.posclusterslabelmat==1);
subplot(3,3,9);ft_singleplotTFR(cfg,statMem);
tit=title(['pval=' num2str(statMem.posclusters(1).prob)]);set(findobj(tit,'type','text'),'FontSize',10);
xlabel('AMY Phase (Hz)','FontSize',10);
ylabel('HC Amplitude (Hz)','FontSize',10);
colormap(gca,'default')
h = colorbar;
h.FontSize = 8;
ylabel(h, 'T-score')
%% Hippocampus memory cluster
load(fullfile(oripath,'data2github','stat_memK_hc.mat'))
mask = sum(squeeze(MemK_hc.mask),2);
fvec = MemK_hc.freq(mask>0);
foi = [min(fvec) max(fvec)];

h1 = Figure(1,'size',[60 30],'fontSizePrint',20);

cfg=[];
cfg.frequency = foi;
cfg.avgoverfreq = 'yes';
cfg.latency = [3 12];
cfg.avgovertime = 'yes';
ef = ft_selectdata(cfg,efreq);
nf = ft_selectdata(cfg,nfreq);

O = teg_repeated_measures_ANOVA([ef.powspctrm,nf.powspctrm], 2, {'emotion'});
titO = ['F_(_' num2str(O.R(1,2)) '_,_' num2str(O.R(1,3)) '_)=' num2str(O.R(1,1)) ', p=' num2str(O.R(1,4))];

x = ones(9,2).*[1:2];
y = [ef.powspctrm nf.powspctrm];

subplot(121);beeswarm(x(:),y(:),'sort_style','up','dot_size',1,'overlay_style','sd','colormap',[0.5430 0 0; 0 0 0.4430]);
a = get(gca,'XTickLabels');
set(gca,'XTickLabel',{'Emotional','Neutral'},'fontsize',8,'FontWeight','bold')
h=get(gca);
ylim([0 1.8e-3]);
ylabel('Modulation Index');title('low\gamma_H_C - \theta_A_M_Y');

load(fullfile(oripath,'data2github','Stat_intK_amy_8.mat'))
mask = sum(squeeze(Int_amy.mask),2);
fvec = Int_amy.freq(mask>0);
foi = [min(fvec) max(fvec)];

cfg=[];
cfg.frequency = foi;
cfg.avgoverfreq = 'yes';
cfg.latency = [3 12];
cfg.avgovertime = 'yes';
ef = ft_selectdata(cfg,efreq);
nf = ft_selectdata(cfg,nfreq);

O = teg_repeated_measures_ANOVA([ef.powspctrm,nf.powspctrm], 2, {'emotion'});
titO = ['F_(_' num2str(O.R(1,2)) '_,_' num2str(O.R(1,3)) '_)=' num2str(O.R(1,1)) ', p=' num2str(O.R(1,4))];

x = ones(9,2).*[1:2];
y = [ef.powspctrm nf.powspctrm];

subplot(122);beeswarm(x(:),y(:),'sort_style','up','dot_size',1,'overlay_style','sd','colormap',[0.5430 0 0; 0 0 0.4430]);
a = get(gca,'XTickLabels');
set(gca,'XTickLabel',{'Emotional','Neutral'},'fontsize',8,'FontWeight','bold')
h=get(gca);
ylim([0 1.8e-3]);
title('high\gamma_H_C - \theta_A_M_Y');

clear eR_amy eKF_amy nR_amy nKF_amy e_amy n_amy r_amy f_amy data* h h1...
      eR_hc eKF_hc nR_hc nKF_hc e_hc n_hc r_hc f_hc
