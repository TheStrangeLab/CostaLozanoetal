
%% This script reproduce the Entorhinal (EC) and Perirhinal (PRh)  time frequency results and statistic. 
% It produces Supplementary Figure 8 b, c, d 

clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20190203\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','TF','ECPR','TimeFreq_cond_ECPRhsubjects.mat'))

%%
subjects=[13 15 34 5 500 600 8 10];
 
%%
 for j= 1:length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
 TFs_Multitaper_baseline_EKF(j) = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_eK(j).values, TFs_Multitaper_baseline_average_eF(j).values)
 TFs_Multitaper_baseline_NKF(j) = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_nK(j).values, TFs_Multitaper_baseline_average_nF(j).values)
 end
 %%
  for j= 1:length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
 TFs_Multitaper_baseline_KF(j) = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_eK(j).values, TFs_Multitaper_baseline_average_eF(j).values,TFs_Multitaper_baseline_average_nK(j).values, TFs_Multitaper_baseline_average_nF(j).values)
 TFs_Multitaper_baseline_K(j) = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_eK(j).values, TFs_Multitaper_baseline_average_nK(j).values)

 end
 %%
for v=1:length(subjects)
cfg = [];
TFeR(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_average_eR(v).values);
TFnR(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_average_nR(v).values);
TFeF(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_average_eF(v).values);
TFnF(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_average_nF(v).values);
TFeK(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_average_eK(v).values);
TFnK(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_average_nK(v).values);
TFeKF(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_EKF(v));
TFnKF(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_NKF(v));
end

%%
for v=1:length(subjects)
cfg = []
TFKF(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_KF(v));
TFK(v) = ft_freqdescriptives(cfg, TFs_Multitaper_baseline_K(v));
end
%% 
for v= 1:length(subjects);
        cfg = [];
        cfg.appenddim ='rpt'
        cfg.parameter = 'powspctrm';
        TFs_Multitaper_bipolar_e(v).values = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_eR(v).values, TFs_Multitaper_baseline_average_eF(v).values);
        TFs_Multitaper_bipolar_n(v).values = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_nR(v).values, TFs_Multitaper_baseline_average_nF(v).values);

end

%%
 for v= 1:length(subjects);
 cfg = [];
 cfg.appenddim ='rpt'
 cfg.parameter = 'powspctrm';
 TFs_Multitaper_bipolar_rem(v).values = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_nR(v).values, TFs_Multitaper_baseline_average_eR(v).values);
 TFs_Multitaper_bipolar_forg(v).values = ft_appendfreq(cfg, TFs_Multitaper_baseline_average_nF(v).values, TFs_Multitaper_baseline_average_eF(v).values);
 end


%%
for v= 1:length(subjects)
cfg = [];
TFs_e(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_e(v).values);
TFs_n(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_n(v).values);
TFs_rem(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_rem(v).values);
TFs_forg(v) = ft_freqdescriptives(cfg,TFs_Multitaper_bipolar_forg(v).values);
end

%%
close all

cfg = [];
cfg.figure ='gcf'
cfg.xlim =[-0.5 1.5];
cfg.ylim =[35 150];

for j= 1:length(subjects)
    h = figure;set(h, 'Position', get(0, 'Screensize'))
    subplot(2,3,1)
        cfg.zlim = [-0.5 0.5];
        cfg.parameter  = 'powspctrm'
        ft_singleplotTFR(cfg,TFeR(j))
        title(['TFs er, Subject ' num2str(subjects(j))]);
    subplot(2,3,2)
        cfg.parameter  = 'powspctrm'
        ft_singleplotTFR(cfg,TFeK(j))
        title(['TFs ek, Subject ' num2str(subjects(j))]);
        hold on;
    subplot(2,3,3)
        ft_singleplotTFR(cfg,TFeF(j))
        title(['TFs ef, Subject ' num2str(subjects(j))]);
        hold on
    subplot(2,3,4)
        cfg.parameter  = 'powspctrm'
        ft_singleplotTFR(cfg,TFnR(j))
        title(['TFs nr, Subject ' num2str(subjects(j))]);
        hold on;
    subplot(2,3,5)
        ft_singleplotTFR(cfg,TFnK(j))
        title(['TFs nk, Subject ' num2str(subjects(j))]);
        hold on
        subplot(2,3,6)
        cfg.parameter  = 'powspctrm'
        ft_singleplotTFR(cfg,TFnF(j))
        title(['TFs nf, Subject ' num2str(subjects(j))]);
     
end
%%

cfg = [];
cfg.channel = {'ParaHcEPC'};
cfg.keepindividual = 'yes';

GTFs_eR = ft_freqgrandaverage(cfg,TFeR(1),TFeR(2), TFeR(3),TFeR(4),TFeR(5),TFeR(6),TFeR(7),TFeR(8));
GTFs_nR = ft_freqgrandaverage(cfg,TFnR(1),TFnR(2), TFnR(3),TFnR(4),TFnR(5),TFnR(6),TFnR(7),TFnR(8));
GTFs_eKF = ft_freqgrandaverage(cfg,TFeKF(1),TFeKF(2), TFeKF(3),TFeKF(4),TFeKF(5),TFeKF(6),TFeKF(7),TFeKF(8));
GTFs_nKF = ft_freqgrandaverage(cfg,TFnKF(1),TFnKF(2), TFnKF(3),TFnKF(4),TFnKF(5),TFnKF(6),TFnKF(7),TFnKF(8));

%%

GTFs_diff_URUKF = GTFs_eR
GTFs_diff_URUKF.powspctrm = GTFs_eR.powspctrm - GTFs_eKF.powspctrm

GTFs_diff_NRNKF = GTFs_nR
GTFs_diff_NRNKF.powspctrm = GTFs_nR.powspctrm - GTFs_nKF.powspctrm

%% Supplementary Fig.8b

cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];
cfg.ylim = [35 150];
cfg.zlim = [0 0.2];
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,GTFs_eR); 
tit=title('eR');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,GTFs_nR);
tit=title('nR');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,3) ; ft_singleplotTFR(cfg,GTFs_eKF); 
tit=title('eKF');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,GTFs_nKF);
tit=title('nKF');set(findobj(tit,'type','text'),'FontSize',36);

%% STAT
cfg = [];
cfg.method           = 'ft_statistics_montecarlo'; 
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusterstatistic = 'maxsum';

t1=0
t2=1
f1=35;
f2=150;

cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

%cfg.minnbchan      = 2;
cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;%
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 'all';%10000;

cfg.neighbours = []; %only one channel -> fieldtrip recognizes time-freq-neighbours
cfg.ivar = 1;
cfg.uvar = 2;

 Nsubjects= size(GTFs_eR.powspctrm,1);
 design = zeros(2,2*Nsubjects);
 
 for i = 1:Nsubjects
     design(2,i) = i;
 end

 for i = 1:Nsubjects
     design(2,Nsubjects+i) = i;
 end
 
design(1,1:Nsubjects)= 1;
design(1,Nsubjects+1:2*Nsubjects) = 2;

cfg.design = design;

IntK = ft_freqstatistics(cfg,GTFs_diff_URUKF, GTFs_diff_NRNKF);

%%
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1]%'maxmin';
cfg.ylim = [35 150]% 'maxmin';
cfg.zlim = [-3 3];

E = IntK;
E.powspctrm = IntK.stat.*IntK.mask;
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('IntK'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16) ;
%% supplementary Fig 8c
h=figure;set(h, 'Position', get(0, 'Screensize'))

logRelative8 =IntK
logRelative8.mask =IntK.mask;
logRelative8.powspctrm =IntK.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative8);

%% find the range
maskt= sum(squeeze(IntK.mask),1);
tvec =IntK.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(IntK.mask),2);
fvec =IntK.freq(maskf>0);
foi=[min(fvec) max(fvec)]
%% mean gamma power
t = toi
f = foi

pt1 = nearest(GTFs_eR.time,t(1));
pt2 = nearest(GTFs_eR.time,t(2));
pf1 = nearest(GTFs_eR.freq,f(1));
pf2 = nearest(GTFs_eR.freq,f(2));

% create matrix of conditions 
mat(:,1) = squeeze(mean(mean(GTFs_eR.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(GTFs_eKF.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,3) = squeeze(mean(mean(GTFs_nR.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,4) = squeeze(mean(mean(GTFs_nKF.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
%%
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'));

x = [ones(8,1) ones(8,1)*2 ones(8,1)*3 ones(8,1)*4];
y = [mat(:,1) mat(:,2) mat(:,3) mat(:,4)];
figure;beeswarm(x(:),y(:),'sort_style','up','dot_size',4,'overlay_style','sd','colormap',[1 0 0; 0 0 0; 0 0 1; 0.5 0.5 0.5])
ylim([-0.3 0.5]);
xticklabels({[],'eR',[],'eKF',[],'nR',[],'nKF'})






