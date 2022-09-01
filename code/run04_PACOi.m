% This script computes the phase to amplitude coupling opposition index (PACOi)
% between amygdala phase and hippocampus amplitude contacts under the eR vs eKF
% contrast. Fig 3a, Supplementary Fig. 19

clear all,close all,clc
rng(144);

mlist={'s6 right', 's13', 's15', 's16', 's16 right', 's25', 's32', 's10', 's33'};%;

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','pac'))
addpath(fullfile(oripath,'Costalozanoetal','code','pacoi'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','circstat-matlab'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','teg_repeated_measures_ANOVA'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','violinPlot'))
addpath(fullfile(oripath,'Costalozanoetal','code','pacoi'))

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

for subj = 1:9
  cfg= [];
  eKF_amy{subj} = ft_appenddata(cfg, eK_amy{subj}, eF_amy{subj});
  nKF_amy{subj} = ft_appenddata(cfg, nK_amy{subj}, nF_amy{subj});
  eKF_hc{subj} = ft_appenddata(cfg, eK_hc{subj}, eF_hc{subj});
  nKF_hc{subj} = ft_appenddata(cfg, nK_hc{subj}, nF_hc{subj});
end
clear eK_amy eF_amy nK_amy nF_amy eK_hc eF_hc nK_hc nF_hc

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


% trials: patients(9,1) x conditions (eR_hc, eKF_hc, nR_hc, nKF_hc)
trl=[];
pval=zeros(9,2);

nperm =1000;
[aheR_x, aheKF_x, ahnR_x, ahnKF_x] = deal(nan(9,size(position,2),size(flow,2),size(fhigh,2)));
[p1, pz1, pos1, pos_a1, pos_b1, p2, pz2, pos2, pos_a2, pos_b2] = deal(nan(9,size(flow,2),size(fhigh,2)));
trl2=[];

for s=[1:9]
  aeR = raw2data(eR_amy{1,s},'rpt_chan_time');
  aeR = squeeze(aeR.trial);
  heR = raw2data(eR_hc{1,s}, 'rpt_chan_time');
  heR = squeeze(heR.trial);
  [mi,meanamp1] = ModIndex_v2_ft(aeR,heR,srate,flow,flow_bw,fhigh,fhigh_bw,position);

  mi2 = permute(mi,[3 1 2]);
  k2sel1 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,1)=size(k2sel1,1);
  MeanAmp = (mean(meanamp1(:,:,:,k2sel1),4));
  eRAm(:,:,:,s) = MeanAmp;
  eRMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];

  meanamp1 = meanamp1(:,:,:,k2sel1);
  siz = size(meanamp1);
  aheR_c = nan(siz(4),siz(3),siz(1),siz(2));
  for f1=1:siz(1)
    for f2=1:siz(2)
      for k=1:siz(4)
        aheR_c(k,:,f1,f2) = squeeze(meanamp1(f1,f2,:,k))'.*exp(sqrt(-1).*position);
      end
    end
  end
  mi=[];

  % amygdala-to-hippocampus emotional forgotten
  aeKF = raw2data(eKF_amy{1,s},'rpt_chan_time');
  aeKF = squeeze(aeKF.trial);
  heKF = raw2data(eKF_hc{1,s}, 'rpt_chan_time');
  heKF = squeeze(heKF.trial);

  [mi,meanamp2] = ModIndex_v2_ft(aeKF,heKF,srate,flow,flow_bw,fhigh,fhigh_bw,position);

  mi2 = permute(mi,[3 1 2]);
  k2sel2 = find(sum(isnan(mi2(:,:)),2)==0);
  trl(s,2)=size(k2sel2,1);

  MeanAmp = (mean(meanamp2(:,:,:,k2sel2),4));
  eFAm(:,:,:,s) = MeanAmp;
  eFMIx(:,:,s)=(log(nbin)-(-sum((MeanAmp./sum(MeanAmp,3)).*log((MeanAmp./sum(MeanAmp,3))),3)))./log(nbin);
  mi=[];

  meanamp2 = meanamp2(:,:,:,k2sel2);
  siz = size(meanamp2);
  aheKF_c = nan(siz(4),siz(3),siz(1),siz(2));
  for f1=1:siz(1)
    for f2=1:siz(2)
      for k=1:siz(4)
        aheKF_c(k,:,f1,f2) = squeeze(meanamp2(f1,f2,:,k))'.*exp(sqrt(-1).*position);
      end
    end
  end
  mi=[];

  aheR_cx = squeeze(mean(aheR_c,2));
  aheKF_cx = squeeze(mean(aheKF_c,2));

  for f1=1:size(flow,2)
    for f2=1:size(fhigh,2)
      [p1(s,f1,f2), pz1(s,f1,f2), pos1(s,f1,f2), pos_a1(s,f1,f2), pos_b1(s,f1,f2)]=...
      ppoi(squeeze(aheR_cx(:,f1,f2)), squeeze(aheKF_cx(:,f1,f2)), nperm);
    end
  end
  [size(aheR_cx,1) size(aheKF_cx,1)]

  mask1 = double(squeeze(pz1(s,:,:))<0.05)';
  mask1(~mask1) = 0.2;

  figure(1);subplot(3,3,s);imagesc(flow, fhigh, squeeze(pos1(s,:,:))',[0 1.2]);
  alpha(mask1);colormap parula
  axis xy;

  % now filter the bar graph
  aheR_x(s,:,:,:) = abs(squeeze(mean(aheR_c,1)));
  aheKF_x(s,:,:,:) = abs(squeeze(mean(aheKF_c,1)));

  figure(2);
  subplot(3,3,s);bar(position,mean(aheR_x(s,:,(squeeze(pz1(s,:,:))<0.05)'),3),'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',2);xlim([position(1) position(end)]);
  subplot(3,3,s);hold all;bar(position,mean(aheKF_x(s,:,(squeeze(pz1(s,:,:))<0.05)'),3),'FaceColor','none','EdgeColor',[0 0.5 1],'LineWidth',2);xlim([position(1) position(end)]);
end

%% plot eR vs eKF
x_dbin = mean([position(1:end-1);position(2:end)],1);
x_dbincorr =[x_dbin 2.9845];
pos_i = [1:2:18];
pos_b = pos_i+1;

eMa=[];
nMa=[];
Figure(7,'size',[120   80],'fontSizeScreen',10);
for s=[1:9]

  mask = double(squeeze(pz1(s,:,:))<0.05)';
  mask(~mask) = 0.3;

  figure(7);
  subplot(3,6,pos_i(s));imagesc(flow, fhigh, squeeze(pos1(s,:,:))',[0 0.3]);colormap('default');
  alpha(mask);title(['Patient nr ' num2str(s)])
  axis xy;
  xlabel('AMY Phase (Hz)','FontSize',8);
  ylabel('HC Amplitude (Hz)','FontSize',8);
  if s==1;
    colorbar;
    h = colorbar;
    h.FontSize = 8;
    ylabel(h, 'PACOi eR vs eKF')
  end

  m1 = sum(squeeze(pz1(s,:,:))<0.05,2)';

%   angbin1b = mean(aheR_x(s,:,(squeeze(pz1(s,:,:))<0.05)'),3);
%   angbin2b = mean(aheKF_x(s,:,(squeeze(pz1(s,:,:))<0.05)'),3);
%------------------------------------------
masksig = (squeeze(pz1(s,:,:))<0.05);
ang1 = squeeze(aheR_x(s,:,:,:));
ang2 = squeeze(aheKF_x(s,:,:,:));

angbin1=[];
angbin2=[];
[r,c]=find(masksig);
for ii=1:size(r,1)
  angbin1(:,ii) = ang1(:,r(ii),c(ii));
  angbin2(:,ii) = ang2(:,r(ii),c(ii));
end
angbin1 = mean(angbin1,2)';
angbin2 = mean(angbin2,2)';
%------------------------------------------
[angbin1a, angbin2a2, flag(s),angbin2a] = realign2phases(angbin1,angbin2,x_dbincorr);
eMa(s,:) = angbin1a;
nMa(s,:) = angbin2a2;

subplot(3,6,pos_b(s));bar(position,zscore(angbin1a,[],2),'FaceColor','none','EdgeColor',[1 0 0],'LineWidth',1);xlim([position(1) position(end)]);hold all;
subplot(3,6,pos_b(s));bar(position,zscore(angbin2a2,[],2),'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1);xlim([position(1) position(end)]);
ylim([-3 3]);
xlabel('AMY phase (rad)','FontSize',8);
ylabel('HC amplitude (z)','FontSize',8);
end

%% compute vector strength and angle
eMc = mean(eMa.*exp(sqrt(-1)*x_dbincorr),2);
nMc = mean(nMa.*exp(sqrt(-1)*x_dbincorr),2);

E_ph = angle(mean(eMc./abs(eMc)));
N_ph = angle(mean(nMc./abs(nMc)));

E_ab = abs(mean(eMc./abs(eMc)));
N_ab = abs(mean(nMc./abs(nMc)));

rE = mean(zscore(eMa,[],2));
rN = mean(zscore(nMa,[],2));

rErose = mean(eMa);
rNrose = mean(nMa);

%% histogram
x = linspace(-0.5,0.5,20);
w=2*pi*1;
W=[w-w/8];
[param1,y_est1,y_resid1,err_rms1] = sinefit(rE,x,W);
[param2,y_est2,y_resid2,err_rms2] = sinefit(rN,x,W);

h1 = Figure(1,'size',[120   120],'fontSizeScreen',10);
subplot(221);bar(x_dbincorr,(rE),'FaceColor','none','EdgeColor',[1 0 0],'LineWidth',1);
subplot(221);hold all;bar(x_dbincorr,(rN),'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1);
subplot(221);plot(E_ph,1.4,'vr','MarkerSize',8,'MarkerFaceColor','r')
subplot(221);plot(N_ph,1.4,'vk','MarkerSize',8,'MarkerFaceColor','k')
xlabel('AMY_\theta phase (rad)');
ylabel('HC_\gamma amplitude (z)');

subplot(221);plot(x_dbincorr,y_est1,'Color',[1 0 0],'LineWidth',2);
subplot(221);hold all;plot(x_dbincorr,y_est2,'Color',[0 0 0],'LineWidth',2);
xlim([-pi-0.1 pi+0.1]);
set(gca, 'box','off')
ll=legend({'eR','eKF','eR angle','eKF angle'},'location','northwest');
set(gca, 'box','off')
set(ll, 'box','off')
%% pacoi roseplot

angle_bins      = deg2rad(-180:18:180);
angle_centers  	= transpose(movmean(angle_bins, 2, 'endpoints', 'discard'));

clear dt
dt.figEx = [2 5 12 8]
dt.boxEx = [4 1.35 4.25 4.25]
dt.visible = 'on'
dt.angLabels{1} = '0°'
dt.angLabels{2} = '90°'
dt.angLabels{3} = '180°'
dt.angLabels{4}= '-90°'
dt.bReverseY = 1
dt.angleCenters = angle_centers
dt.angularRes = 20
dt.FR = rE';
dt.FR1 =  rN';

h= plot_rosegroup_pacoi(dt)
set(h, 'PaperPositionMode', 'auto');
%print(h, 'groupres_pacoi_erekf', '-dtiff', '-r450');


