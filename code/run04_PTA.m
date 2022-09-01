% This script computes the peak trigger averages (PTA) within and between regions.
% Fig 2ghi, Supplementary Fig. 15,16,17
clear all,close all,clc
rng(144);

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20190203\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','pac'))
addpath(fullfile(oripath,'Costalozanoetal','code','pacoi'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','circstat-matlab'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','teg_repeated_measures_ANOVA'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','violinPlot'))
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
%% peak triggered averages for Emotion vs Neutral:
% Hippocampus gamma peaks - Amygdala raw traces

srate = 500;
dt = 1/srate;

win = 60;
minpeakdistance = 100;

fhigh = [50 75];
fhigh_bw = 0;

cfg=[];
cfg.bpfilter = 'yes';
cfg.bpfreq = [fhigh(1)-(fhigh_bw/2) fhigh(2)+(fhigh_bw/2)];
cfg.bpfilttype = 'fir';
cfg.bpfiltdir = 'twopass';
cfg.demean = 'yes';

avg_haE = [];
avg_haN = [];
for subj = 1:9
  gamma_ehc = ft_preprocessing(cfg,data_emo_hc{1,subj});
  gamma_nhc = ft_preprocessing(cfg,data_neu_hc{1,subj});

  ehc_g = raw2data(gamma_ehc,'rpt_chan_time');
  ehc_g = squeeze(ehc_g.trial);
  nhc_g = raw2data(gamma_nhc,'rpt_chan_time');
  nhc_g = squeeze(nhc_g.trial);

  eamy_r = raw2data(data_emo_amy{1,subj},'rpt_chan_time');
  eamy_r = squeeze(eamy_r.trial);
  namy_r = raw2data(data_neu_amy{1,subj},'rpt_chan_time');
  namy_r = squeeze(namy_r.trial);

  avg1 = [];
  avg2 = [];
  loc1 = [];
  loc2 = [];
  cnt1 = [];
  cnt2 = [];
  for trl = 1:size(eamy_r,1)
    [g,lag,l] = peaktriggeredavg(eamy_r(trl,:),ehc_g(trl,:),win,minpeakdistance);
    avg1 = cat(1,g,avg1);
    loc1 = cat(2,l,loc1);
    cnt1 = [cnt1 size(l,2)];
    clear g lag l
  end
  for trl = 1:size(namy_r,1)
    [g,lag,l] = peaktriggeredavg(namy_r(trl,:),nhc_g(trl,:),win,minpeakdistance);
    avg2 = cat(1,g,avg2);
    loc2 = cat(2,l,loc2);
    cnt2 = [cnt1 size(l,2)];
    clear g lag l
  end

  avg_haE(subj,:) = detrend(mean(zscore(avg1,[],2))')';
  loc_ha{subj} = loc1;
  cnt_ha{subj} = cnt1;

  avg_haN(subj,:) = detrend(mean(zscore(avg2,[],2))')';
  loc_hh{subj} = loc2;
  cnt_hh{subj} = cnt2;
end

rem   = [1 0 0];
kforg = [0 0 0];
emo   = [0.7 0.0 0.0];
neu   = [0   0   0.4];

pad = 512;
normalize = 'sum01';

%close all
x = (-win:win)./srate;
lat = nearest(x,[-0.09 0.09]);
h2=Figure(2,'size',[100 100],'fontSizeScreen',10);
for s=[1:9]

  z = zscore(avg_haE(s,:));
  [val,idx]=max(abs(z(1,lat(1):lat(2))));

  [p1,f1]=fftpsd(avg_haE(s,:), dt, size(x,2), pad, normalize);
  [p2,f2]=fftpsd(avg_haN(s,:), dt, size(x,2), pad, normalize);
  lag(s)= x(idx+lat(1));

  [vl,id]=max(p1);
  fpeak =  f2(id);

  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  [param1,y_est1,y_resid1,err_rms1] = sinefit(avg_haE(s,:),x,W,1e-4,64);

  [vl,id]=max(p2);
  fpeak =  f2(id);

  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  [param2,y_est2,y_resid2,err_rms2] = sinefit(avg_haN(s,:),x,W,1e-4,64);

  figure(2);sh2 = subplot(3,3,s);pos=sh2.Position;sh2.Position(4)=pos(4)*0.5;

  [xr(s,:),l] = xcorr(avg_haE(s,:),avg_haN(s,:),'coeff');
  [p3,f3]=fftpsd(xr(s,:), dt, size(xr,2), pad, normalize);
  psd(s,:) = p3;

  flim = nearest(f3,[0 15]);
  [vl,id]=max(p3(flim(1):flim(2)));
  p3 = p3./max(p3);
  fpeak = f3(id);
  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  lim = (1/(fpeak))/2;
  xlm = nearest(l/500,[-lim lim]);
  xl0 = nearest(l/500,0);

  [param3,y_est3,y_resid3,err_rms3] = sinefit(xr(s,:),l/500,W,1e-4,64);

  [rmaxnone,imaxnone]=max(xr(s,xlm(1):xlm(2)));
  tau1(s) = l(xlm(1)+imaxnone)/500;

  figure(2);subplot(sh2);hold all;plot([0 0],[-0.6 0.6],'--k','linewidth',1);

  x1=l/500;
  y=xr(s,:);
  z=y;

  long_dash = xr(s,:);
  long_dash([1:20 41:60 81:100 121:140 161:180 201:220])=NaN;

  figure(2);subplot(sh2);hold all;plot([0 0],[-0.9 0.9],'--k','linewidth',1);
  figure(2);subplot(sh2);hold all;ax=plot(l/500,xr(s,:),'linewidth',1.35,'Color',emo);
  figure(2);subplot(sh2);hold all;ax=plot(l/500,long_dash,'linewidth',1.35,'Color',neu);
  figure(2);subplot(sh2);hold all;plot(l/500,y_est3,'Color',[0 0 0],'linewidth',2);
  %xlabel('lag (s)');
  %ylabel('corr (r)');
  xlim([-0.2 0.2]);ylim([-0.7 0.7]);
  ax = gca;
  ax.XAxisLocation = 'origin';
  ax.YAxisLocation = 'origin';
  set(gca, 'box' , 'off');
  
  npos = sh2.Position;
  npos(4) = npos(4)*0.3;
  npos(2) = npos(2)-(npos(4)+0.025);
  nax =axes('Position',npos);subplot(nax);box(nax,'off');
  subplot(nax);hold all;ax=plot(f3,p3,'linewidth',1,'Color',[0 0 0]);
  xlim([0 25]);ylim([0 1.2]);
  nax =axes('Position',npos);subplot(nax);box(nax,'off');
  subplot(nax);hold all;ax=plot(f3,p3,'linewidth',1,'Color',[0 0 0]);
  xlim([0 25]);ylim([0 1.2]);

end

%% ------------------------------------------------------------------------
% peak triggered averages for Emotion vs Neutral:
% Hippocampus gamma peaks - Hippocampus raw traces
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [fhigh(1)-(fhigh_bw/2) fhigh(2)+(fhigh_bw/2)];
cfg.bpfilttype = 'fir';
cfg.bpfiltdir = 'twopass';
cfg.demean = 'yes';

avg_haE = [];
avg_haN = [];

for subj = 1:9
  gamma_ehc = ft_preprocessing(cfg,data_emo_hc{1,subj});
  gamma_nhc = ft_preprocessing(cfg,data_neu_hc{1,subj});

  ehc_g = raw2data(gamma_ehc,'rpt_chan_time');
  ehc_g = squeeze(ehc_g.trial);
  nhc_g = raw2data(gamma_nhc,'rpt_chan_time');
  nhc_g = squeeze(nhc_g.trial);

  ehc_r = raw2data(data_emo_hc{1,subj},'rpt_chan_time');
  ehc_r = squeeze(ehc_r.trial);
  nhc_r = raw2data(data_neu_hc{1,subj},'rpt_chan_time');
  nhc_r = squeeze(nhc_r.trial);

  avg1 = [];
  avg2 = [];
  loc1 = [];
  loc2 = [];
  cnt1 = [];
  cnt2 = [];
  for trl = 1:size(ehc_r,1)
    [g,lag,l] = peaktriggeredavg(ehc_r(trl,:),ehc_g(trl,:),win,minpeakdistance);
    avg1 = cat(1,g,avg1);
    loc1 = cat(2,l,loc1);
    cnt1 = [cnt1 size(l,2)];
    clear g lag l
  end
  for trl = 1:size(nhc_r,1)
    [g,lag,l] = peaktriggeredavg(nhc_r(trl,:),nhc_g(trl,:),win,minpeakdistance);
    avg2 = cat(1,g,avg2);
    loc2 = cat(2,l,loc2);
    cnt2 = [cnt1 size(l,2)];
    clear g lag l
  end

  avg_haE(subj,:) = detrend(mean(zscore(avg1,[],2))')';
  loc_ha{subj} = loc1;
  cnt_ha{subj} = cnt1;

  avg_haN(subj,:) = detrend(mean(zscore(avg2,[],2))')';
  loc_hh{subj} = loc2;
  cnt_hh{subj} = cnt2;
end

%%
%close all
h3=Figure(3,'size',[100 100],'fontSizeScreen',10);
x = (-win:win)./srate;
lat = nearest(x,[-0.09 0.09]);
xr=[];
for s=[1:9]

  z = avg_haE(s,:);
  [val,idx]=max(abs(z(1,lat(1):lat(2))));

  [p1,f1]=fftpsd(avg_haE(s,:), 1/srate, size(x,2), pad, normalize);
  [p2,f2]=fftpsd((avg_haN(s,:)), 1/srate, size(x,2), pad, normalize);
  lag(s)= x(idx+lat(1));

  [vl,id]=max(p1);
  fpeak =  f2(id);

  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  [param1,y_est1,y_resid1,err_rms1] = sinefit(avg_haE(s,:),x,W,1e-4,64);

  flim = nearest(f1,[0 15]);
  [vl,id]=max(p2(flim(1):flim(2)));
  fpeak =  f2(id);

  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  [param2,y_est2,y_resid2,err_rms2] = sinefit((avg_haN(s,:)),x,W,1e-4,64);

  figure(3);sh3 = subplot(3,3,s);pos=sh3.Position;sh3.Position(4)=pos(4)*0.5;
  [xr(s,:),l] = xcorr(avg_haE(s,:),(avg_haN(s,:)),'coeff');
  [p3,f3]=fftpsd(xr(s,:), 1/srate, size(xr,2), pad, normalize);
  psd(s,:) = p3;

  flim = nearest(f3,[0 15]);
  [vl,id]=max(p3(flim(1):flim(2)));
  p3 = p3./max(p3);
  fpeak = f3(id);
  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  lim = (1/(fpeak))/2;
  xlm = nearest(l/500,[-lim lim]);
  xl0 = nearest(l/500,0);

  [param3,y_est3,y_resid3,err_rms3] = sinefit(xr(s,:),l/500,W,1e-4,64);

  [rmaxnone,imaxnone]=max(xr(s,xlm(1):xlm(2)));
  tau1(s) = l(xlm(1)+imaxnone)/500;

  long_dash = xr(s,:);
  long_dash([1:20 41:60 81:100 121:140 161:180 201:220])=NaN;

  figure(3);subplot(sh3);hold all;plot([0 0],[-0.9 0.9],'--k','linewidth',1);
  figure(3);subplot(sh3);hold all;ax=plot(l/500,xr(s,:),'linewidth',1.35,'Color',emo);
  figure(3);subplot(sh3);hold all;ax=plot(l/500,long_dash,'linewidth',1.35,'Color',neu);
  figure(3);subplot(sh3);hold all;plot(l/500,y_est3,'Color',[0 0 0],'linewidth',2);
  xlabel('lag (s)');
  ylabel('corr (r)');
  xlim([-0.2 0.2]);ylim([-0.78 0.78]);
  ax = gca;
  ax.XAxisLocation = 'origin';
  ax.YAxisLocation = 'origin';
  set(gca, 'box' , 'off');

  npos = sh3.Position;
  npos(4) = npos(4)*0.3;
  npos(2) = npos(2)-(npos(4)+0.025);
  nax =axes('Position',npos);subplot(nax);box(nax,'off');
  subplot(nax);hold all;ax=plot(f3,p3,'linewidth',1,'Color',[0 0 0]);
  xlim([0 25]);ylim([0 1.2]);
  nax =axes('Position',npos);subplot(nax);box(nax,'off');
  subplot(nax);hold all;ax=plot(f3,p3,'linewidth',1,'Color',[0 0 0]);
  xlim([0 25]);ylim([0 1.2]);
  hold on
end

%% peak triggered averages for eR vs eF:
% Hippocampus gamma peaks - Amygdala raw traces
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [fhigh(1)-(fhigh_bw/2) fhigh(2)+(fhigh_bw/2)];
cfg.bpfilttype = 'fir';
cfg.bpfiltdir = 'twopass';
cfg.demean = 'yes';

avg_haR = [];
avg_haF = [];
for subj = 1:9
  gamma_eRhc = ft_preprocessing(cfg,eR_hc{1,subj});
  gamma_eFhc = ft_preprocessing(cfg,eKF_hc{1,subj});

  eRhc_g = raw2data(gamma_eRhc,'rpt_chan_time');
  eRhc_g = squeeze(eRhc_g.trial);
  eFhc_g = raw2data(gamma_eFhc,'rpt_chan_time');
  eFhc_g = squeeze(eFhc_g.trial);

  eRamy_r = raw2data(eR_amy{1,subj},'rpt_chan_time');
  eRamy_r = squeeze(eRamy_r.trial);
  eFamy_r = raw2data(eKF_amy{1,subj},'rpt_chan_time');
  eFamy_r = squeeze(eFamy_r.trial);

  avg1 = [];
  avg2 = [];
  loc1 = [];
  loc2 = [];
  cnt1 = [];
  cnt2 = [];
  for trl = 1:size(eRhc_g,1)
    [g,lag,l] = peaktriggeredavg(eRamy_r(trl,:),eRhc_g(trl,:),win,minpeakdistance);
    avg1 = cat(1,g,avg1);
    loc1 = cat(2,l,loc1);
    cnt1 = [cnt1 size(l,2)];
    clear g lag l
  end
  for trl = 1:size(eFhc_g,1)
    [g,lag,l] = peaktriggeredavg(eFamy_r(trl,:),eFhc_g(trl,:),win,minpeakdistance);
    avg2 = cat(1,g,avg2);
    loc2 = cat(2,l,loc2);
    cnt2 = [cnt1 size(l,2)];
    clear g lag l
  end

  avg_haR(subj,:) = detrend(mean(zscore(avg1,[],2))')';
  loc_ha{subj} = loc1;
  cnt_ha{subj} = cnt1;

  avg_haF(subj,:) = detrend(mean(zscore(avg2,[],2))')';
  loc_hh{subj} = loc2;
  cnt_hh{subj} = cnt2;
end


%close all
h4=Figure(4,'size',[100 100],'fontSizeScreen',10);
x = (-win:win)./srate;
 lat = nearest(x,[-0.09 0.09]);
xr=[];
for s=[1:9]

  z = (avg_haR(s,:));
  [val,idx]=max(abs(z(1,lat(1):lat(2))));

  [p1,f1]=fftpsd((avg_haR(s,:)), 1/srate, size(x,2), pad, normalize);
  [p2,f2]=fftpsd((avg_haF(s,:)), 1/srate, size(x,2), pad, normalize);
  lag(s)= x(idx+lat(1));


 if s==14 || s==19;
  flim = nearest(f1,[0 10]);
    [vl,id]=max(p1(flim(1):flim(2)));
    fpeak =  f2(id);
  else
    [vl,id]=max(p1);
    fpeak =  f2(id);
  end
  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  [param1,y_est1,y_resid1,err_rms1] = sinefit((avg_haR(s,:)),x,W,1e-4,64);

  if s== 2|| s==19;
  flim = nearest(f1,[0 10]);
    [vl,id]=max(p2(flim(1):flim(2)));
    fpeak =  f2(id);
  else
    [vl,id]=max(p2);
    fpeak =  f2(id);
  end

  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  [param2,y_est2,y_resid2,err_rms2] = sinefit((avg_haF(s,:)),x,W,1e-4,64);

  figure(4);sh4 = subplot(3,3,s);pos=sh4.Position;sh4.Position(4)=pos(4)*0.5;
  
  [xr(s,:),l] = xcorr((avg_haR(s,:)),(avg_haF(s,:)),'coeff');
  [p3,f3]=fftpsd(xr(s,:), 1/srate, size(xr,2), pad, normalize);
  psd(s,:) = p3;

  flim = nearest(f3,[0 15]);
  [vl,id]=max(p3(flim(1):flim(2)));
  p3 = p3./max(p3);
  fpeak = f3(id);
  w=2*pi*fpeak;
  W=[w-w/4,w+w/4];
  lim = (1/(fpeak))/2;
  xlm = nearest(l/500,[-lim lim]);
  xl0 = nearest(l/500,0);

  [param3,y_est3,y_resid3,err_rms3] = sinefit(xr(s,:),l/500,W,1e-4,64);

  [rmaxnone,imaxnone]=max(xr(s,xlm(1):xlm(2)));
  tau1(s) = l(xlm(1)+imaxnone)/500;

  long_dash = xr(s,:);
  long_dash([1:20 41:60 81:100 121:140 161:180 201:220])=NaN;

  figure(4);subplot(sh4);hold all;ax=plot(l/500,xr(s,:),'linewidth',1.35,'Color',rem);
  figure(4);subplot(sh4);hold all;ax=plot(l/500,long_dash,'linewidth',1.35,'Color',kforg);
  figure(4);subplot(sh4);hold all;plot(l/500,y_est3,'Color',[0 0 0],'linewidth',2);
  %xlabel('lag (s)');
  %ylabel('corr (r)');
  xlim([-0.2 0.2]);ylim([-0.7 0.7]);
  ax = gca;
  ax.XAxisLocation = 'origin';
  ax.YAxisLocation = 'origin';
  set(gca, 'box' , 'off');

  npos = sh4.Position;
  npos(4) = npos(4)*0.3;
  npos(2) = npos(2)-(npos(4)+0.025);
  nax =axes('Position',npos);subplot(nax);box(nax,'off');
  subplot(nax);hold all;ax=plot(f3,p3,'linewidth',1,'Color',[0 0 0]);
  xlim([0 25]);ylim([0 1.2]);
  nax =axes('Position',npos);subplot(nax);box(nax,'off');
  subplot(nax);hold all;ax=plot(f3,p3,'linewidth',1,'Color',[0 0 0]);
  xlim([0 25]);ylim([0 1.2]);
end
