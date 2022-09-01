%% This script reproduce the hippocampal time frequency results and statistic. 
% It also produces Figure 1 g,h,i and Supplementary Fig.5 and 7

clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20190203\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','TF','Hippo','TimeFreq_cond_9sj.mat'))

subjects=[60 13 15 16 160 25 32 10 33] % 10 is z1 right
%% Do Avaragae with baseline as average of channels

TF_eR = [];
for j = 1:length(subjects)
   
        TF_eR(j).values = eR(j).values;
        if size(eR(j).values.powspctrm,2) ~= 1 %channels
            TF_eR(j).values.powspctrm = mean(eR(j).values.powspctrm,2); 
        end
        TF_eR(j).values.label = {'Hippo'};
        size(TF_eR(j).values.powspctrm)
end


TF_nR = [];
for j = 1:length(subjects)
   
        TF_nR(j).values = nR(j).values;
        if size(nR(j).values.powspctrm,2) ~= 1 
            TF_nR(j).values.powspctrm = mean(nR(j).values.powspctrm,2); 
        end
        TF_nR(j).values.label = {'Hippo'};
        size(TF_nR(j).values.powspctrm)
    end


TF_eF = [];
for j = 1:length(subjects)
   
        TF_eF(j).values = eF(j).values;
        if size(eF(j).values.powspctrm,2) ~= 1 
            TF_eF(j).values.powspctrm = mean(eF(j).values.powspctrm,2); 
        end
        TF_eF(j).values.label = {'Hippo'};
        size(TF_eF(j).values.powspctrm)
    end


TF_nF = [];
for j = 1:length(subjects)
   
        TF_nF(j).values = nF(j).values;
        if size(nF(j).values.powspctrm,2) ~= 1 
            TF_nF(j).values.powspctrm = mean(nF(j).values.powspctrm,2); 
        end
        TF_nF(j).values.label = {'Hippo'};
        size(TF_nF(j).values.powspctrm)
    end


TF_eK = [];
for j = 1:length(subjects)
   
        TF_eK(j).values = eK(j).values;
        if size(eK(j).values.powspctrm,2) ~= 1 
            TF_eK(j).values.powspctrm = mean(eK(j).values.powspctrm,2); 
        end
        TF_eK(j).values.label = {'Hippo'};
        size(TF_eK(j).values.powspctrm)
    end


TF_nK = [];
for j = 1:length(subjects)
   
        TF_nK(j).values = nK(j).values;
        if size(nK(j).values.powspctrm,2) ~= 1 
            TF_nK(j).values.powspctrm = mean(nK(j).values.powspctrm,2); 
        end
        TF_nK(j).values.label = {'Hippo'};
        size(TF_nK(j).values.powspctrm)
    end


TF_eKF = [];
for j = 1:length(subjects)
   
        TF_eKF(j).values = eKF(j).values;
        if size(eKF(j).values.powspctrm,2) ~= 1 
            TF_eKF(j).values.powspctrm = mean(eKF(j).values.powspctrm,2); 
        end
        TF_eKF(j).values.label = {'Hippo'};
        size(TF_eKF(j).values.powspctrm)
    end


TF_nKF = [];
for j = 1:length(subjects)
   
        TF_nKF(j).values = nKF(j).values;
        if size(nKF(j).values.powspctrm,2) ~= 1 
            TF_nKF(j).values.powspctrm = mean(nKF(j).values.powspctrm,2);  
        end
        TF_nKF(j).values.label = {'Hippo'};
        size(TF_nKF(j).values.powspctrm)
    end

%% main effects

 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        TFs_Multitaper_bipolar_Unpl(j).values = ft_appendfreq(cfg, TF_eR(j).values, TF_eF(j).values)
        TFs_Multitaper_bipolar_UnplKF(j).values = ft_appendfreq(cfg, TF_eR(j).values, TF_eKF(j).values)
        TFs_Multitaper_bipolar_Neu(j).values = ft_appendfreq(cfg, TF_nR(j).values, TF_nF(j).values)
        TFs_Multitaper_bipolar_NeuKF(j).values = ft_appendfreq(cfg, TF_nR(j).values, TF_nKF(j).values)
 end
 
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        TFs_Multitaper_bipolar_R(j).values = ft_appendfreq(cfg, TF_eR(j).values, TF_nR(j).values);
        TFs_Multitaper_bipolar_F(j).values = ft_appendfreq(cfg, TF_eF(j).values, TF_nF(j).values);
        TFs_Multitaper_bipolar_K(j).values = ft_appendfreq(cfg, TF_eK(j).values, TF_nK(j).values);
        TFs_Multitaper_bipolar_KF(j).values = ft_appendfreq(cfg, TF_eKF(j).values, TF_nKF(j).values);
 end
 

%%
 for v= 1:length(subjects)
cfg = [];
TF_avg_chann_trl_eR(v) = ft_freqdescriptives(cfg, TF_eR(v).values);
TF_avg_chann_trl_eK(v) = ft_freqdescriptives(cfg, TF_eK(v).values);
TF_avg_chann_trl_eF(v) = ft_freqdescriptives(cfg, TF_eF(v).values);
TF_avg_chann_trl_eKF(v) = ft_freqdescriptives(cfg, TF_eKF(v).values);

TF_avg_chann_trl_nR(v) = ft_freqdescriptives(cfg, TF_nR(v).values);
TF_avg_chann_trl_nK(v) = ft_freqdescriptives(cfg, TF_nK(v).values);
TF_avg_chann_trl_nF(v) = ft_freqdescriptives(cfg, TF_nF(v).values);
TF_avg_chann_trl_nKF(v) = ft_freqdescriptives(cfg, TF_nKF(v).values);
end


%%
for v= 1:length(subjects)
cfg = [];
TFs_Multitaper_unpl(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_Unpl(v).values);
TFs_Multitaper_unplkf(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_UnplKF(v).values);
TFs_Multitaper_neu(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_Neu(v).values);
TFs_Multitaper_neukf(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_NeuKF(v).values);
end


%%
for v= 1:length(subjects)
cfg = [];
TFs_Multitaper_rem(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_R(v).values);
TFs_Multitaper_forg(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_F(v).values);
TFs_Multitaper_kf(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_KF(v).values);
TFs_Multitaper_k(v) = ft_freqdescriptives(cfg, TFs_Multitaper_bipolar_K(v).values);
end

%% Compute grand average 

cfg = [];
cfg.channel = {'Hippo'};
cfg.keepindividual = 'yes';

GTFs_Multitaper_baseline_eR = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eR(1),TF_avg_chann_trl_eR(2),TF_avg_chann_trl_eR(3),TF_avg_chann_trl_eR(4), TF_avg_chann_trl_eR(5),TF_avg_chann_trl_eR(6),TF_avg_chann_trl_eR(7),TF_avg_chann_trl_eR(8),TF_avg_chann_trl_eR(9));
GTFs_Multitaper_baseline_nR = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nR(1),TF_avg_chann_trl_nR(2), TF_avg_chann_trl_nR(3),TF_avg_chann_trl_nR(4),TF_avg_chann_trl_nR(5),TF_avg_chann_trl_nR(6),TF_avg_chann_trl_nR(7),TF_avg_chann_trl_nR(8),TF_avg_chann_trl_nR(9));
GTFs_Multitaper_baseline_eF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eF(1),TF_avg_chann_trl_eF(2),TF_avg_chann_trl_eF(3),TF_avg_chann_trl_eF(4), TF_avg_chann_trl_eF(5),TF_avg_chann_trl_eF(6),TF_avg_chann_trl_eF(7),TF_avg_chann_trl_eF(8),TF_avg_chann_trl_eF(9));
GTFs_Multitaper_baseline_nF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nF(1),TF_avg_chann_trl_nF(2),TF_avg_chann_trl_nF(3),TF_avg_chann_trl_nF(4), TF_avg_chann_trl_nF(5),TF_avg_chann_trl_nF(6),TF_avg_chann_trl_nF(7),TF_avg_chann_trl_nF(8),TF_avg_chann_trl_nF(9));
GTFs_Multitaper_baseline_eK = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eK(1),TF_avg_chann_trl_eK(2),TF_avg_chann_trl_eK(3),TF_avg_chann_trl_eK(4), TF_avg_chann_trl_eK(5),TF_avg_chann_trl_eK(6),TF_avg_chann_trl_eK(7),TF_avg_chann_trl_eK(8),TF_avg_chann_trl_eK(9));
GTFs_Multitaper_baseline_nK = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nK(1),TF_avg_chann_trl_nK(2),TF_avg_chann_trl_nK(3),TF_avg_chann_trl_nK(4), TF_avg_chann_trl_nK(5),TF_avg_chann_trl_nK(6),TF_avg_chann_trl_nK(7),TF_avg_chann_trl_nK(8),TF_avg_chann_trl_nK(9));
GTFs_Multitaper_baseline_eKF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eKF(1),TF_avg_chann_trl_eKF(2),TF_avg_chann_trl_eKF(3),TF_avg_chann_trl_eKF(4), TF_avg_chann_trl_eKF(5),TF_avg_chann_trl_eKF(6),TF_avg_chann_trl_eKF(7),TF_avg_chann_trl_eKF(8),TF_avg_chann_trl_eKF(9));
GTFs_Multitaper_baseline_nKF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nKF(1),TF_avg_chann_trl_nKF(2),TF_avg_chann_trl_nKF(3),TF_avg_chann_trl_nKF(4), TF_avg_chann_trl_nKF(5),TF_avg_chann_trl_nKF(6),TF_avg_chann_trl_nKF(7),TF_avg_chann_trl_nKF(8),TF_avg_chann_trl_nKF(9));

%% calculate grandaverage unpl and neu
cfg = [];
cfg.channel = {'Hippo'};
cfg.keepindividual = 'yes';

GTFs_Multitaper_baseline_Unpl = ft_freqgrandaverage(cfg,TFs_Multitaper_unpl(1),TFs_Multitaper_unpl(2),TFs_Multitaper_unpl(3),TFs_Multitaper_unpl(4),TFs_Multitaper_unpl(5),TFs_Multitaper_unpl(6), TFs_Multitaper_unpl(7),TFs_Multitaper_unpl(8),TFs_Multitaper_unpl(9));%
GTFs_Multitaper_baseline_Unplkf = ft_freqgrandaverage(cfg,TFs_Multitaper_unplkf(1),TFs_Multitaper_unplkf(2),TFs_Multitaper_unplkf(3),TFs_Multitaper_unplkf(4),TFs_Multitaper_unplkf(5),TFs_Multitaper_unplkf(6), TFs_Multitaper_unplkf(7),TFs_Multitaper_unplkf(8),TFs_Multitaper_unplkf(9));%
GTFs_Multitaper_baseline_Neu = ft_freqgrandaverage(cfg,TFs_Multitaper_neu(1),TFs_Multitaper_neu(2),TFs_Multitaper_neu(3),TFs_Multitaper_neu(4),TFs_Multitaper_neu(5),TFs_Multitaper_neu(6), TFs_Multitaper_neu(7),TFs_Multitaper_neu(8),TFs_Multitaper_neu(9));%
GTFs_Multitaper_baseline_Neukf = ft_freqgrandaverage(cfg,TFs_Multitaper_neukf(1),TFs_Multitaper_neukf(2),TFs_Multitaper_neukf(3),TFs_Multitaper_neukf(4),TFs_Multitaper_neukf(5),TFs_Multitaper_neukf(6), TFs_Multitaper_neukf(7),TFs_Multitaper_neukf(8),TFs_Multitaper_neukf(9));%

%% calculate grandaverage R Forg
cfg = [];
cfg.channel = {'Hippo'};
cfg.keepindividual = 'yes';

GTFs_Multitaper_baseline_R = ft_freqgrandaverage(cfg,TFs_Multitaper_rem(1),TFs_Multitaper_rem(2),TFs_Multitaper_rem(3),TFs_Multitaper_rem(4),TFs_Multitaper_rem(5),TFs_Multitaper_rem(6), TFs_Multitaper_rem(7),TFs_Multitaper_rem(8),TFs_Multitaper_rem(9));
GTFs_Multitaper_baseline_F = ft_freqgrandaverage(cfg,TFs_Multitaper_forg(1),TFs_Multitaper_forg(2),TFs_Multitaper_forg(3),TFs_Multitaper_forg(4),TFs_Multitaper_forg(5),TFs_Multitaper_forg(6), TFs_Multitaper_forg(7),TFs_Multitaper_forg(8),TFs_Multitaper_forg(9));
GTFs_Multitaper_baseline_K = ft_freqgrandaverage(cfg,TFs_Multitaper_k(1),TFs_Multitaper_k(2),TFs_Multitaper_k(3),TFs_Multitaper_k(4),TFs_Multitaper_k(5),TFs_Multitaper_k(6), TFs_Multitaper_k(7),TFs_Multitaper_k(8),TFs_Multitaper_k(9));
GTFs_Multitaper_baseline_KF = ft_freqgrandaverage(cfg,TFs_Multitaper_kf(1),TFs_Multitaper_kf(2),TFs_Multitaper_kf(3),TFs_Multitaper_kf(4),TFs_Multitaper_kf(5),TFs_Multitaper_kf(6), TFs_Multitaper_kf(7),TFs_Multitaper_kf(8),TFs_Multitaper_kf(9));

%% 
GTFs_diff_URUF = GTFs_Multitaper_baseline_eR
GTFs_diff_URUF.powspctrm = GTFs_Multitaper_baseline_eR.powspctrm - GTFs_Multitaper_baseline_eF.powspctrm

GTFs_diff_NRNF = GTFs_Multitaper_baseline_nR
GTFs_diff_NRNF.powspctrm = GTFs_Multitaper_baseline_nR.powspctrm - GTFs_Multitaper_baseline_nF.powspctrm

%% 
GTFs_diff_URUKF = GTFs_Multitaper_baseline_eR
GTFs_diff_URUKF.powspctrm = GTFs_Multitaper_baseline_eR.powspctrm - GTFs_Multitaper_baseline_eKF.powspctrm

GTFs_diff_NRNKF = GTFs_Multitaper_baseline_nR
GTFs_diff_NRNKF.powspctrm = GTFs_Multitaper_baseline_nR.powspctrm - GTFs_Multitaper_baseline_nKF.powspctrm

%% Plot Grand Average Multitaper WITH baseline High frequency
% Fig 1g
cfg = []; 
cfg.xlim = [-.5 1.5];%
cfg.ylim = [35 150];%'maxmin';
cfg.zlim = [0 0.2];%[0 0.2];%
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_eR);
tit=title('eR');set(findobj(tit,'type','text'),'FontSize',36); 
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_nR); 
tit=title('nR');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,3) ; ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_eKF);
tit=title('eKF');set(findobj(tit,'type','text'),'FontSize',36);   
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_nKF); 
tit=title('nKF');set(findobj(tit,'type','text'),'FontSize',36);
 
%% Plot Grand Average Multitaper WITH baseline High frequency
%Supplementary Fig.5b
cfg = []; 
cfg.xlim = [-.5 1.5];
cfg.ylim = [0 34];
cfg.zlim = [0 1];
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_eR); 
tit=title('eR');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_nR);
tit=title('nR');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,3) ; ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_eKF); 
tit=title('eKF');set(findobj(tit,'type','text'),'FontSize',36);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_nKF);
tit=title('nKF');set(findobj(tit,'type','text'),'FontSize',36);

%% STAT

cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=35;
f2=150;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

cfg.tail             = 0; 
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 10000;

cfg.neighbours = [];
cfg.ivar = 1;
cfg.uvar = 2;

cfg.design = [ones(1,size(subjects,2)) ones(1,size(subjects,2)).*2;[1:size(subjects,2)] [1:size(subjects,2)]];

MemK_hc = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_R, GTFs_Multitaper_baseline_KF);
EmoK_hc = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_Unplkf, GTFs_Multitaper_baseline_Neukf)
IntK_hc = ft_freqstatistics(cfg,GTFs_diff_URUKF, GTFs_diff_NRNKF);
%% plot the significant effect
oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','TF','Hippo','stat_memK_hc.mat'))
% find the range
maskt= sum(squeeze(MemK_hc.mask),1);
tvec = MemK_hc.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(MemK_hc.mask),2);
fvec = MemK_hc.freq(maskf>0);
foi=[min(fvec) max(fvec)]
%% mean gamma

t= toi
f = foi

pt1 = nearest(GTFs_diff_URUKF.time,t(1));
pt2 = nearest(GTFs_diff_URUKF.time,t(2));
pf1 = nearest(GTFs_diff_URUKF.freq,f(1));
pf2 = nearest(GTFs_diff_URUKF.freq,f(2));

% create matrix of conditions 
mat(:,1) = squeeze(mean(mean(GTFs_Multitaper_baseline_eR.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(GTFs_Multitaper_baseline_eKF.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,3) = squeeze(mean(mean(GTFs_Multitaper_baseline_nR.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,4) = squeeze(mean(mean(GTFs_Multitaper_baseline_nKF.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));

%%
 [He,Pe,CIe,STATSe] = ttest(mat(:,1),mat(:,2));
 [Hn,Pn,CIn,STATSn] = ttest(mat(:,3),mat(:,4));

 %% calculate d cohen
 der = mean(mat(:,1))
 dekf = mean(mat(:,2))
 dnr = mean(mat(:,3))
 dnkf = mean(mat(:,4))
 
 appmat_emo = [mat(:,1)-mat(:,2)];%[mat(:,1);mat(:,2)]
 stdemo = std(appmat_emo)
 
 appmat_neu = [mat(:,3)-mat(:,4)];%[mat(:,3);mat(:,4)]
 stdneu = std(appmat_neu)
 
 de = (der-dekf)/stdemo
 dn= (dnr-dnkf)/stdneu
%% Fig 1i
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'));

x = [ones(9,1) ones(9,1)*2 ones(9,1)*3 ones(9,1)*4 ones(9,1)*5 ones(9,1)*6];
y = [mat(:,1) nan(9,1) mat(:,2:3) nan(9,1) mat(:,4)];
figure;beeswarm(x(:),y(:),'sort_style','up','dot_size',4,'overlay_style','sd','colormap',[1 0 0; 1 1 1;0 0 0; 0 0 1; 1 1 1;0.5 0.5 0.5])
ylim([-0.2 0.5]);
xticklabels({'eR',[],'eKF','nR',[],'nKF'})

grid off

%%
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1.5]
cfg.ylim = [35 150]
cfg.zlim = [-3 3];
cfg.maskstyle= 'outline'

M = MemK_hc;
M.powspctrm = MemK_hc.stat.*MemK_hc.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('Main effect of memory K'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16) ; 

%% Fig 1h
h=figure;set(h, 'Position', get(0, 'Screensize'))

logRelative2 = MemK_hc
logRelative2.mask = MemK_hc.mask;
logRelative2.powspctrm = MemK_hc.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative2);
%%
t= toi
f = foi

pt1 = nearest(GTFs_diff_URUKF.time,t(1));
pt2 = nearest(GTFs_diff_URUKF.time,t(2));
pf1 = nearest(GTFs_diff_URUKF.freq,f(1));
pf2 = nearest(GTFs_diff_URUKF.freq,f(2));

% create matrix of conditions 
matme(:,1) = squeeze(mean(mean(GTFs_Multitaper_baseline_R.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
matme(:,2) = squeeze(mean(mean(GTFs_Multitaper_baseline_KF.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
%% Supplementary Fig. 7a

cfg = []; 
cfg.xlim = [-.5 1.5];
cfg.ylim = [35 150];
cfg.zlim = [0 0.2];
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,3,1) ; ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_eR); 
tit=title('eR');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,3,4); ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_nR);
tit=title('nR');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,3,2) ; ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_eK);
tit=title('eK');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,3,5); ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_nK);
tit=title('nK');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,3,3) ; ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_eF);  
tit=title('eF');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,3,6); ft_singleplotTFR(cfg,GTFs_Multitaper_baseline_nF); 
tit=title('nF');set(findobj(tit,'type','text'),'FontSize',12);

%% STAT

cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=35;
f2=150;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

cfg.tail             = 0; 
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;%
cfg.numrandomization = 10000;

cfg.neighbours = []; 
cfg.ivar = 1;
cfg.uvar = 2;

cfg.design = [ones(1,9) ones(1,9).*2;[1:9] [1:9]];

Mem_hcRvsF = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_R, GTFs_Multitaper_baseline_F);
Mem_hcRvsK = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_R, GTFs_Multitaper_baseline_K);

%%
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1.5]
cfg.ylim = [35 150]
cfg.zlim = [-3 3];
cfg.maskstyle= 'outline'

M = Mem_hcRvsF;
Mem_hcRvsF.mask = Mem_hcRvsF.posclusterslabelmat ==1
M.powspctrm = Mem_hcRvsF.stat.*Mem_hcRvsF.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('Main effect of memory R vs F'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16) ; 

%% Supplementary Fig.7c
h=figure;set(h, 'Position', get(0, 'Screensize'))

logRelative2 = Mem_hcRvsF
logRelative2.mask = Mem_hcRvsF.mask;
logRelative2.powspctrm = Mem_hcRvsF.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative2);

%%
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1.5]
cfg.ylim = [35 150]
cfg.zlim = [-3 3];
cfg.maskstyle= 'outline'

M = Mem_hcRvsK;
M.powspctrm = Mem_hcRvsK.stat.*Mem_hcRvsK.mask;
M.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,M); colorbar;  tit=title('Main effect of memory R vs K'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16) ; 

%% Supplementary Fig. 7b
h=figure;set(h, 'Position', get(0, 'Screensize'))

logRelative2 = Mem_hcRvsK
logRelative2.mask = Mem_hcRvsK.mask;
logRelative2.powspctrm = Mem_hcRvsK.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative2);