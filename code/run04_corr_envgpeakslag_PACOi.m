% This script computes the correlation between PACOi and the amplitude
% envelope cross-correlation method. Fig 5b, Supplementary Fig. 20

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
addpath(fullfile(oripath,'Costalozanoetal','code','pacoi'))

load(fullfile(oripath,'data2github','DataCleanAH_RKF.mat'))
% column meaning: eR_hc eKF_hc nR_hc nKF_hc
load(fullfile(oripath,'data2github','trials_patients.mat'))

toi = [0.41 1.1];

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
  data_emo_amy{subj} = ft_appenddata(cfg, eR_amy{subj}, eK_amy{subj}, eF_amy{subj});
  data_neu_amy{subj} = ft_appenddata(cfg, nR_amy{subj}, nK_amy{subj}, nF_amy{subj});
  data_emo_hc{subj}  = ft_appenddata(cfg, eR_hc{subj},  eK_hc{subj},  eF_hc{subj});
  data_neu_hc{subj}  = ft_appenddata(cfg, nR_hc{subj},  nK_hc{subj},  nF_hc{subj});
end

for subj = 1:9
  cfg= [];
  eKF_amy{subj} = ft_appenddata(cfg, eK_amy{subj}, eF_amy{subj});
  nKF_amy{subj} = ft_appenddata(cfg, nK_amy{subj}, nF_amy{subj});
  eKF_hc{subj} = ft_appenddata(cfg, eK_hc{subj}, eF_hc{subj});
  nKF_hc{subj} = ft_appenddata(cfg, nK_hc{subj}, nF_hc{subj});
end
clear eK_amy eF_amy nK_amy nF_amy eK_hc eF_hc nK_hc nF_hc

%% peak triggered averages for eR vs eF:
% Hippocampus gamma peaks - Amygdala raw traces
win = 15;
minpeakdistance = 50;

cfg=[];
cfg.bpfilter = 'yes';
cfg.bpfreq = [60 120];
cfg.bpfilttype = 'fir';
cfg.bpfiltdir = 'twopass';
cfg.demean = 'yes';

avg_haR = [];
avg_haF = [];
sem_haR = [];
sem_haF = [];
ntrials = [];
Cer=[];
Cef=[];
Figure(7,'size',[120   80],'fontSizeScreen',8);
for subj = 1:9
  gamma_eRhc = ft_preprocessing(cfg,eR_hc{1,subj});
  gamma_eFhc = ft_preprocessing(cfg,eKF_hc{1,subj});

  eRhc_g = raw2data(gamma_eRhc,'rpt_chan_time');
  eRhc_g = squeeze(eRhc_g.trial);
  eFhc_g = raw2data(gamma_eFhc,'rpt_chan_time');
  eFhc_g = squeeze(eFhc_g.trial);

  eRhc_r = raw2data(eR_hc{1,subj},'rpt_chan_time');
  eRhc_r = squeeze(eRhc_r.trial);
  eFhc_r = raw2data(eKF_hc{1,subj},'rpt_chan_time');
  eFhc_r = squeeze(eFhc_r.trial);


  gamma_eRamyg = ft_preprocessing(cfg,eR_amy{1,subj});
  gamma_eFamyg = ft_preprocessing(cfg,eKF_amy{1,subj});
  eRamy_g = raw2data(gamma_eRamyg,'rpt_chan_time');
  eRamy_g = squeeze(eRamy_g.trial);
  eFamy_g = raw2data(gamma_eFamyg,'rpt_chan_time');
  eFamy_g = squeeze(eFamy_g.trial);

  eRamy_r = raw2data(eR_amy{1,subj},'rpt_chan_time');
  eRamy_r = squeeze(eRamy_r.trial);
  eFamy_r = raw2data(eKF_amy{1,subj},'rpt_chan_time');
  eFamy_r = squeeze(eFamy_r.trial);

  [rg1,rg2,rr,rxl,lag,lr] = peaktriggeredxcorravg(eRhc_g,eRamy_g,eRhc_g,win,minpeakdistance);
  [fg1,fg2,fr,fxl,lag,lf] = peaktriggeredxcorravg(eFhc_g,eFamy_g,eFhc_g,win,minpeakdistance);

  err = ft_preproc_hilbert(rr,'abs');
  efr = ft_preproc_hilbert(fr,'abs');

  if subj==8
    % this part of the code was written to select a representative trials
    % to show the amplitude envelope cross-correlation method
    trk = 2;
    npk = 2;
    tlong = eR_amy{1,subj}.time{1};

    figure;
    subplot(2,2,[1 2]);plot(tlong,eRhc_g(trk,:) ,'Color',[0 0 .8]);xlim([tlong(1) tlong(end)]);hold all;ylim([-10 10])
    subplot(2,2,[1 2]);plot(tlong,eRamy_g(trk,:),'Color',[0 0.5 0]);xlim([tlong(1) tlong(end)]);hold all;ylim([-10 10])
    subplot(2,2,[1 2]);plot(tlong(lr{trk}),ones(size(lr{trk})).*8,'.k');
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');

    for p=1:size(lr{trk},2)
      lm = lag+lr{trk}(p);
      hold all; subplot(2,2,[1 2]);plot(tlong(lag+lr{trk}(p)),ones(size(lm)).*-5,'-k');
    end
    box off;

    sit = cumsum(cellfun(@length, lr));
    tshort = tlong(lag+lr{trk}(npk));
    subplot(223);plot(tshort,[rg1(sit(trk-1)+npk,:)'],'color',[0 0 .8],'linewidth',1);xlim([tshort(1) tshort(end)]);hold all;ylim([-10 10])
    subplot(223);plot(tshort,[rg2(sit(trk-1)+npk,:)'],'color',[0 0.5 0],'linewidth',1);xlim([tshort(1) tshort(end)]);ylim([-10 10])
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    box off;

    xrtime = rxl./500;
    subplot(224);plot(xrtime,[rr(sit(trk-1)+npk,:)'],'color',[0.5 0.5 0.5],'linewidth',1);xlim([xrtime(1) xrtime(end)]);hold all;
    h=subplot(224);plot(xrtime,[err(sit(trk-1)+npk,:)'],'color',[1 0 0],'linewidth',2);xlim([xrtime(1) xrtime(end)])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Lag (s)');
    ylabel('R_x_y');
    set (h, 'box' , 'off');

    trk = 4;
    npk = 1;
    figure;
    subplot(2,2,[1 2]);plot(tlong,eFhc_g(trk,:) ,'Color',[0 0 .8]);xlim([tlong(1) tlong(end)]);hold all;ylim([-10 10])
    subplot(2,2,[1 2]);plot(tlong,eFamy_g(trk,:),'Color',[0 0.5 0]);xlim([tlong(1) tlong(end)]);hold all;ylim([-10 10])
    subplot(2,2,[1 2]);plot(tlong(lf{trk}),ones(size(lf{trk})).*8,'.k');
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    box off;
    for p=1:size(lf{trk},2)
      lm = lag+lf{trk}(p);
      hold all; subplot(2,2,[1 2]);plot(tlong(lag+lf{trk}(p)),ones(size(lm)).*-5,'-k');
    end

    sit = cumsum(cellfun(@length, lf));
    tshort = tlong(lag+lf{trk}(npk));
    subplot(223);plot(tshort,[fg1(sit(trk-1)+npk,:)'],'color',[0 0 .8],'linewidth',1);xlim([tshort(1) tshort(end)]);hold all;ylim([-10 10])
    subplot(223);plot(tshort,[fg2(sit(trk-1)+npk,:)'],'color',[0 0.5 0],'linewidth',1);xlim([tshort(1) tshort(end)]);ylim([-10 10])
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    box off;
    xrtime = rxl./500;
    subplot(224);plot(xrtime,[fr(sit(trk-1)+npk,:)'],'color',[0.5 0.5 0.5],'linewidth',1);xlim([xrtime(1) xrtime(end)]);hold all;
    h=subplot(224);plot(xrtime,[efr(sit(trk-1)+npk,:)'],'color',[0 0 0],'linewidth',2);xlim([xrtime(1) xrtime(end)])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Lag (s)');
    ylabel('R_x_y');
    set (h, 'box' , 'off');
    
    p = mean(err,2);
    r_idx = p<prctile(p,92) & p>prctile(p,60);

    p = mean(efr,2);
    f_idx = p<prctile(p,92) & p>prctile(p,60);

    figure;
    subplot(221);plot(xrtime,[err(r_idx,:)'],'color',[1 0 0 0.2],'linewidth',1);xlim([xrtime(1) xrtime(end)]);hold all;
    h=subplot(221);plot(xrtime,mean(err(r_idx,:))','color',[1 0 0],'linewidth',2);xlim([xrtime(1) xrtime(end)])
    ylim([-0.3 5.5]);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Lag (s)');
    ylabel('env(R_x_y)');
    set (h, 'box' , 'off');

    subplot(222);plot(xrtime,[efr(f_idx,:)'],'color',[0 0 0 0.2],'linewidth',1);xlim([xrtime(1) xrtime(end)]);hold all;
    h=subplot(222);plot(xrtime,[mean(efr(f_idx,:))'],'color',[0 0 0],'linewidth',2);xlim([xrtime(1) xrtime(end)])
    ylim([-0.3 5.5]);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Lag (s)');
    ylabel('env(R_x_y)');
    set (h, 'box' , 'off');

    long_dash = mean(err(r_idx,:))-mean(efr(f_idx,:));
    long_dash([1:10 21:30 41:49])=NaN;

    subplot(223);plot(xrtime,mean(err(r_idx,:))-mean(efr(f_idx,:)),'color',[1 0 0],'linewidth',2);xlim([xrtime(1) xrtime(end)]);hold all;
    h=subplot(223);plot(xrtime,long_dash,'color',[0 0 0],'linewidth',2);xlim([xrtime(1) xrtime(end)]);hold all;
    ylim([-0.3 5.5]);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Lag (s)');
    ylabel('eR - eKF');
    set (h, 'box' , 'off');
  end

  % stanard error of the eman
  d1m = mean(err);    % d1 mean
  d1n = size(err,1);  % d1 num of samples
  d1var = var(err);   % d1 variance
  d1dof = d1n - 1;   % d1 dof

  d2m = mean(efr);    % d2 mean
  d2n = size(efr,1);  % d2 num of samples
  d2var = var(efr);   % d2 variance
  d2dof = d2n - 1;   % d2 dof

  dof = d1n + d2n - 2;  % totoal dof

  % standard error of mean
  dstd = sqrt((d1var*d1dof+d2var*d2dof)/dof);
  dstdm = dstd * sqrt(1/d1n+1/d2n);

  avg_haR(subj,:) = mean(err);
  sem_haR(subj,:) = sem(err);

  avg_haF(subj,:) = mean(efr);
  sem_haF(subj,:) = sem(efr);

  rxl = rxl./500;
  lag = lag./500;
  long_dash = d1m-d2m;
  long_dash([1:10 21:30 41:49])=NaN;

  mm = max(abs(minmax(d1m-d2m+dstdm)));

  figure(7);h=subplot(3,3,subj);fill([rxl rxl(end:-1:1) rxl(1)],[d1m-d2m-dstdm d1m(1,end:-1:1)-d2m(1,end:-1:1)+dstdm(1,end:-1:1)  d1m(:,1)-d2m(:,1)-dstdm(:,1)],[0.5 0.5 0.5],'LineStyle','none');
  xlim([rxl(1) rxl(end)]);
  ylim([-mm(1)-0.2 mm(1)+0.2]);alpha(0.2);hold all;
  figure(7);h=subplot(3,3,subj);plot(rxl,d1m-d2m,'linewidth',1.5,'Color','r');hold all;
  figure(7);h=subplot(3,3,subj);plot(rxl,long_dash,'linewidth',1.5,'Color','k');

  xlabel('Time (s)');
  ylabel('Env. xcorr. (a.u.)')
  title(['Patient ' mlist{subj}],'Fontsize',10)
  ax = gca;
  ax.XAxisLocation = 'origin';
  ax.YAxisLocation = 'origin';
  set (h, 'box' , 'off');
end

%%

t1=[];
t1.label={'chan1'};
t1.time = rxl;
t1.individual(:,1,:) = avg_haR;
t1.dimord = 'subj_chan_time';

t2 = t1;
t2.individual(:,1,:) = avg_haF;

cfg = [];
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusterstatistic = 'maxsum';
cfg.latency          = [];
cfg.avgovertime      = 'no';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;
cfg.neighbours = [];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,9) ones(1,9).*2; [1:9] [1:9]];
stat = ft_timelockstatistics(cfg,t1,t2);
h_out=(stat.posclusterslabelmat==1);

cohensd = mean(avg_haR-avg_haF) ./ std(avg_haR-avg_haF,[],1);

mu = mean(avg_haR-avg_haF);
sd = sem(avg_haR-avg_haF);
long_dash = mu;
long_dash([1:10 21:30 41:49])=NaN;


Figure(2,'size',[80   35],'fontSizeScreen',8);

half_left = nearest(rxl,[rxl(1) 0]);
h_l = rxl(half_left(1):half_left(2));

half_right = nearest(rxl,[0 rxl(end)]);
h_r = rxl(half_right(1):half_right(2));

subplot(121);fill([h_l h_l(end:-1:1) h_l(1)],[mu(half_left(1):half_left(2))-sd(half_left(1):half_left(2)) mu(1,half_left(2):-1:half_left(1))+sd(1,half_left(2):-1:half_left(1))  mu(:,1)-sd(:,1)],[0.5 0.5 0.5],'LineStyle','none');hold all;
subplot(121);fill([h_r h_r(end:-1:1) h_r(1)],[mu(half_right(1):half_right(2))-sd(half_right(1):half_right(2)) mu(1,half_right(2):-1:half_right(1))+sd(1,half_right(2):-1:half_right(1))  mu(:,half_right(1))-sd(:,half_right(1))],[0 1 1],'LineStyle','none');
xlim([rxl(1) rxl(end)]);alpha(0.2);hold all;

h1(1)=subplot(121);plot(rxl,mu,'linewidth',2,'Color','r');hold all;
h1(2)=subplot(121);plot(rxl,long_dash,'linewidth',2,'Color','k')
xlabel('Time (s)');
ylabel('Env. xcorr (a.u.)')
title('xcorr_(_H_i_p_p_o_c_a_m_p_u_s_-_A_m_y_g_d_a_l_a_)','FontSize',10)
subplot(121);hold all;plot(rxl(h_out),h_out(h_out).*0.01,'-b','linewidth',3);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca, 'box' , 'off');
ll=legend(h1([1]),'eR-eKF')
set(ll,'box','off');

dif = avg_haR-avg_haF;
[val,idx2]=max(abs(dif),[],2);


pacoi = load(fullfile(oripath,'data2github','tmp_pacoi_stat9_0.31-1.2_v3.mat'),'pz1','aheR_x','aheKF_x','flow');
delay=[];
nbin = 20; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin+1
  position(j) = -pi+(j-1)*winsize;
end
x_dbin = mean([position(1:end-1);position(2:end)],1);


for s=1:9
  mask = (squeeze(pacoi.pz1(s,:,:))<0.05);
  ang1 = squeeze(pacoi.aheR_x(s,:,:,:));
  ang2 = squeeze(pacoi.aheKF_x(s,:,:,:));

  angbin1=[];
  angbin2=[];
  [r,c]=find(mask);
  for ii=1:size(r,1)
    angbin1(:,ii) = ang1(:,r(ii),c(ii));
    angbin2(:,ii) = ang2(:,r(ii),c(ii));
  end
  angbin1 = mean(angbin1,2)';
  angbin2 = mean(angbin2,2)';

  [angbin1a, angbin2a2,flag] = realign2phases(angbin1,angbin2,x_dbin);

  eMc = mean(angbin1.*exp(sqrt(-1)*x_dbin),2);
  nMc = mean(angbin2.*exp(sqrt(-1)*x_dbin),2);

  po = angle(eMc./abs(eMc)) - angle(nMc./abs(nMc));
  m1 = sum(squeeze(pacoi.pz1(s,:,:))<0.05,2)';

  frrs (s,:) = [sum(m1.*pacoi.flow)./sum(m1) mean(pacoi.flow(m1>0))];
  fr = sum(m1.*pacoi.flow)./sum(m1);
  %  fr = mean(pacoi.flow(m1>0));
  delay(s,:) = po./(2*pi)*(1/fr);
end

type='spearman';
[r,p]=corr(rxl(idx2)',(delay),'type',type);
h1=subplot(122);plot(rxl(idx2)',delay,'or','Markersize',5,'LineWidth',1,'MarkerFaceColor', 'k')
set(h1, 'box' , 'off');
xlabel('xcorr lag (secs)');
ylabel('PACOi delay (secs)');

pbaspect([2 2 2])
hl=lsline;
hl.Color=[0 0 0];
title(['rho =' num2str(r) ', p=' num2str(p)])
xlim([-0.09 0.09]);
ylim([-0.09 0.09]);

nReps = 10000;
perm = zeros(nReps,1);
for i=1:nReps
  [null,index] = sort(rand(size(delay)));
  perm(i) = corr(rxl(idx2)',delay(index),'type',type);
end
p = sum(perm>r)/nReps;

figure;hist(perm,50);
ylim = get(gca,'YLim');
hold all;plot(r*[1,1],ylim,'r-','LineWidth',2);
xlabel("Spearman's \rho")
ylabel('counts')
