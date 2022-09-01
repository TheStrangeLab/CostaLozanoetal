% This script reproduce the amygdala time frequency results and statistic
% reported at Figure 1 d,e,f;

clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20190203\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','TF','Amy','TimeFreq_cond_allamysubjects.mat'))
%% n=17
subjects=[2 4 6 60 13 15 16 160 25 27 32 320 21 210 10 33 34]
mlist={'s2', 's4', 's6', 's60', 's13', 's15', 's16', 's160', 's25', 's27','s32','s320','s21', 's210','sz10','s33','s34'};
behavlist = {'Patient2','Patient4','Patient6','Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient27','Patient32','Patient320','Patient21','Patient210','PatientZ_10','Patient33','Patient34'};

%% Do Avaragae with baseline 
%collapsing channels All
TF_eR = [];
for j = 1:length(subjects)
   
        TF_eR(j).values = eR(j).values;
        if size(eR(j).values.powspctrm,2) ~= 1 %channels
            TF_eR(j).values.powspctrm = mean(eR(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_eR(j).values.label = {'Amy'};
        size(TF_eR(j).values.powspctrm)
end

%% Do Average WITH Baseline correction for Neutral R

%collapsing channels All
TF_nR = [];
for j = 1:length(subjects)
   
        TF_nR(j).values = nR(j).values;
        if size(nR(j).values.powspctrm,2) ~= 1 %channels
            TF_nR(j).values.powspctrm = mean(nR(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_nR(j).values.label = {'Amy'};
        size(TF_nR(j).values.powspctrm)
    end

%% Do Average WITH Baseline correction for ef

%collapsing channels All
TF_eF = [];
for j = 1:length(subjects)
   
        TF_eF(j).values = eF(j).values;
        if size(eF(j).values.powspctrm,2) ~= 1 %channels
            TF_eF(j).values.powspctrm = mean(eF(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_eF(j).values.label = {'Amy'};
        size(TF_eF(j).values.powspctrm)
    end


%% Do Average WITH Baseline correction for nF
%collapsing channels All
TF_nF = [];
for j = 1:length(subjects)
   
        TF_nF(j).values = nF(j).values;
        if size(nF(j).values.powspctrm,2) ~= 1 %channels
            TF_nF(j).values.powspctrm = mean(nF(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_nF(j).values.label = {'Amy'};
        size(TF_nF(j).values.powspctrm)
    end

%% Do Average WITH Baseline correction for ek

%collapsing channels All
TF_eK = [];
for j = 1:length(subjects)
   
        TF_eK(j).values = eK(j).values;
        if size(eK(j).values.powspctrm,2) ~= 1 %channels
            TF_eK(j).values.powspctrm = mean(eK(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_eK(j).values.label = {'Amy'};
        size(TF_eK(j).values.powspctrm)
    end


%% Do Average WITH Baseline correction for nK
%collapsing channels All
TF_nK = [];
for j = 1:length(subjects)
   
        TF_nK(j).values = nK(j).values;
        if size(nK(j).values.powspctrm,2) ~= 1 %channels
            TF_nK(j).values.powspctrm = mean(nK(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_nK(j).values.label = {'Amy'};
        size(TF_nK(j).values.powspctrm)
    end

%% Do Average WITH Baseline correction for ekf

%collapsing channels All
TF_eKF = [];
for j = 1:length(subjects)
   
        TF_eKF(j).values = eKF(j).values;
        if size(eKF(j).values.powspctrm,2) ~= 1 %channels
            TF_eKF(j).values.powspctrm = mean(eKF(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_eKF(j).values.label = {'Amy'};
        size(TF_eKF(j).values.powspctrm)
    end


%% Do Average WITH Baseline correction for Neutral Miss
%collapsing channels All
TF_nKF = [];
for j = 1:length(subjects)
   
        TF_nKF(j).values = nKF(j).values;
        if size(nKF(j).values.powspctrm,2) ~= 1 %channels
            TF_nKF(j).values.powspctrm = mean(nKF(j).values.powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TF_nKF(j).values.label = {'Amy'};
        size(TF_nKF(j).values.powspctrm)
    end

%%

 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        TFs_Multitaper_bipolar_Unpl(j).values = ft_appendfreq(cfg, TF_eR(j).values, TF_eF(j).values)
        TFs_Multitaper_bipolar_UnplKF(j).values = ft_appendfreq(cfg, TF_eR(j).values, TF_eK(j).values, TF_eF(j).values)
        TFs_Multitaper_bipolar_Neu(j).values = ft_appendfreq(cfg, TF_nR(j).values, TF_nF(j).values)
        TFs_Multitaper_bipolar_NeuKF(j).values = ft_appendfreq(cfg, TF_nR(j).values, TF_nK(j).values, TF_nF(j).values)
 end
 
%%

 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        TFs_Multitaper_bipolar_R(j).values = ft_appendfreq(cfg, TF_eR(j).values, TF_nR(j).values);
        TFs_Multitaper_bipolar_F(j).values = ft_appendfreq(cfg, TF_eF(j).values, TF_nF(j).values);
        TFs_Multitaper_bipolar_K(j).values = ft_appendfreq(cfg, TF_eK(j).values, TF_nK(j).values);
        TFs_Multitaper_bipolar_KF(j).values = ft_appendfreq(cfg, TF_eK(j).values, TF_eF(j).values,TF_nK(j).values, TF_nF(j).values);
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

%%
cfg = [];
cfg.channel = {'Amy'};
cfg.keepindividual = 'yes';

Ga_eR = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eR(1),TF_avg_chann_trl_eR(2),TF_avg_chann_trl_eR(3),TF_avg_chann_trl_eR(4),TF_avg_chann_trl_eR(5),TF_avg_chann_trl_eR(6),TF_avg_chann_trl_eR(7),TF_avg_chann_trl_eR(8),TF_avg_chann_trl_eR(9),TF_avg_chann_trl_eR(10),TF_avg_chann_trl_eR(11),TF_avg_chann_trl_eR(12),TF_avg_chann_trl_eR(13),TF_avg_chann_trl_eR(14),TF_avg_chann_trl_eR(16),TF_avg_chann_trl_eR(17));%
Ga_eK = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eK(1),TF_avg_chann_trl_eK(2),TF_avg_chann_trl_eK(3),TF_avg_chann_trl_eK(4),TF_avg_chann_trl_eK(5),TF_avg_chann_trl_eK(6),TF_avg_chann_trl_eK(7),TF_avg_chann_trl_eK(8),TF_avg_chann_trl_eK(9),TF_avg_chann_trl_eK(10),TF_avg_chann_trl_eK(11),TF_avg_chann_trl_eK(12),TF_avg_chann_trl_eK(13),TF_avg_chann_trl_eK(14),TF_avg_chann_trl_eK(16),TF_avg_chann_trl_eK(17));%
Ga_eF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eF(1),TF_avg_chann_trl_eF(2),TF_avg_chann_trl_eF(3),TF_avg_chann_trl_eF(4),TF_avg_chann_trl_eF(5),TF_avg_chann_trl_eF(6),TF_avg_chann_trl_eF(7),TF_avg_chann_trl_eF(8),TF_avg_chann_trl_eF(9),TF_avg_chann_trl_eF(10),TF_avg_chann_trl_eF(11),TF_avg_chann_trl_eF(12),TF_avg_chann_trl_eF(13),TF_avg_chann_trl_eF(14),TF_avg_chann_trl_eF(16),TF_avg_chann_trl_eF(17));%
Ga_nR = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nR(1),TF_avg_chann_trl_nR(2),TF_avg_chann_trl_nR(3),TF_avg_chann_trl_nR(4),TF_avg_chann_trl_nR(5),TF_avg_chann_trl_nR(6),TF_avg_chann_trl_nR(7),TF_avg_chann_trl_nR(8),TF_avg_chann_trl_nR(9),TF_avg_chann_trl_nR(10),TF_avg_chann_trl_nR(11),TF_avg_chann_trl_nR(12),TF_avg_chann_trl_nR(13),TF_avg_chann_trl_nR(14),TF_avg_chann_trl_nR(16),TF_avg_chann_trl_nR(17));%
Ga_nK = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nK(1),TF_avg_chann_trl_nK(2),TF_avg_chann_trl_nK(3),TF_avg_chann_trl_nK(4),TF_avg_chann_trl_nK(5),TF_avg_chann_trl_nK(6),TF_avg_chann_trl_nK(7),TF_avg_chann_trl_nK(8),TF_avg_chann_trl_nK(9),TF_avg_chann_trl_nK(10),TF_avg_chann_trl_nK(11),TF_avg_chann_trl_nK(12),TF_avg_chann_trl_nK(13),TF_avg_chann_trl_nK(14),TF_avg_chann_trl_nK(16),TF_avg_chann_trl_nK(17));%
Ga_nF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nF(1),TF_avg_chann_trl_nF(2),TF_avg_chann_trl_nF(3),TF_avg_chann_trl_nF(4),TF_avg_chann_trl_nF(5),TF_avg_chann_trl_nF(6),TF_avg_chann_trl_nF(7),TF_avg_chann_trl_nF(8),TF_avg_chann_trl_nF(9),TF_avg_chann_trl_nF(10),TF_avg_chann_trl_nF(11),TF_avg_chann_trl_nF(12),TF_avg_chann_trl_nF(13),TF_avg_chann_trl_nF(14),TF_avg_chann_trl_nF(16),TF_avg_chann_trl_nF(17));%
Ga_eKF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_eKF(1),TF_avg_chann_trl_eKF(2),TF_avg_chann_trl_eKF(3),TF_avg_chann_trl_eKF(4),TF_avg_chann_trl_eKF(5),TF_avg_chann_trl_eKF(6),TF_avg_chann_trl_eKF(7),TF_avg_chann_trl_eKF(8),TF_avg_chann_trl_eKF(9),TF_avg_chann_trl_eKF(10),TF_avg_chann_trl_eKF(11),TF_avg_chann_trl_eKF(12),TF_avg_chann_trl_eKF(13),TF_avg_chann_trl_eKF(14),TF_avg_chann_trl_eKF(16),TF_avg_chann_trl_eKF(17));%
Ga_nKF = ft_freqgrandaverage(cfg,TF_avg_chann_trl_nKF(1),TF_avg_chann_trl_nKF(2),TF_avg_chann_trl_nKF(3),TF_avg_chann_trl_nKF(4),TF_avg_chann_trl_nKF(5),TF_avg_chann_trl_nKF(6),TF_avg_chann_trl_nKF(7),TF_avg_chann_trl_nKF(8),TF_avg_chann_trl_nKF(9),TF_avg_chann_trl_nKF(10),TF_avg_chann_trl_nKF(11),TF_avg_chann_trl_nKF(12),TF_avg_chann_trl_nKF(13),TF_avg_chann_trl_nKF(14),TF_avg_chann_trl_nKF(16),TF_avg_chann_trl_nKF(17));%

%load(fullfile(oripath,'data2github','TF','Amy','Ga_conds_16sj.mat'))
%% calculate grandaverage unpl and neu
cfg = [];
cfg.channel = {'Amy'};
cfg.keepindividual = 'yes';

GTFs_Multitaper_baseline_Unpl = ft_freqgrandaverage(cfg,TFs_Multitaper_unpl(1),TFs_Multitaper_unpl(2),TFs_Multitaper_unpl(3),TFs_Multitaper_unpl(4),TFs_Multitaper_unpl(5),TFs_Multitaper_unpl(6), TFs_Multitaper_unpl(7),TFs_Multitaper_unpl(8),TFs_Multitaper_unpl(9),TFs_Multitaper_unpl(10),TFs_Multitaper_unpl(11),TFs_Multitaper_unpl(12),TFs_Multitaper_unpl(13),TFs_Multitaper_unpl(14),TFs_Multitaper_unpl(15), TFs_Multitaper_unpl(16),TFs_Multitaper_unpl(17));
GTFs_Multitaper_baseline_Unplkf = ft_freqgrandaverage(cfg,TFs_Multitaper_unplkf(1),TFs_Multitaper_unplkf(2),TFs_Multitaper_unplkf(3),TFs_Multitaper_unplkf(4),TFs_Multitaper_unplkf(5),TFs_Multitaper_unplkf(6), TFs_Multitaper_unplkf(7),TFs_Multitaper_unplkf(8),TFs_Multitaper_unplkf(9),TFs_Multitaper_unplkf(10),TFs_Multitaper_unplkf(11),TFs_Multitaper_unplkf(12),TFs_Multitaper_unplkf(13),TFs_Multitaper_unplkf(14),TFs_Multitaper_unplkf(15), TFs_Multitaper_unplkf(16),TFs_Multitaper_unplkf(17));
GTFs_Multitaper_baseline_Neu = ft_freqgrandaverage(cfg,TFs_Multitaper_neu(1),TFs_Multitaper_neu(2),TFs_Multitaper_neu(3),TFs_Multitaper_neu(4),TFs_Multitaper_neu(5),TFs_Multitaper_neu(6), TFs_Multitaper_neu(7),TFs_Multitaper_neu(8),TFs_Multitaper_neu(9),TFs_Multitaper_neu(10),TFs_Multitaper_neu(11),TFs_Multitaper_neu(12),TFs_Multitaper_neu(13),TFs_Multitaper_neu(14),TFs_Multitaper_neu(15), TFs_Multitaper_neu(16),TFs_Multitaper_neu(17));
GTFs_Multitaper_baseline_Neukf = ft_freqgrandaverage(cfg,TFs_Multitaper_neukf(1),TFs_Multitaper_neukf(2),TFs_Multitaper_neukf(3),TFs_Multitaper_neukf(4),TFs_Multitaper_neukf(5),TFs_Multitaper_neukf(6), TFs_Multitaper_neukf(7),TFs_Multitaper_neukf(8),TFs_Multitaper_neukf(9),TFs_Multitaper_neukf(10),TFs_Multitaper_neukf(11),TFs_Multitaper_neukf(12),TFs_Multitaper_neukf(13),TFs_Multitaper_neukf(14),TFs_Multitaper_neukf(15), TFs_Multitaper_neukf(16),TFs_Multitaper_neukf(17));

%% calculate grandaverage R Forg
cfg = [];
cfg.channel = {'Amy'};
cfg.keepindividual = 'yes';

GTFs_Multitaper_baseline_R = ft_freqgrandaverage(cfg,TFs_Multitaper_rem(1),TFs_Multitaper_rem(2),TFs_Multitaper_rem(3),TFs_Multitaper_rem(4),TFs_Multitaper_rem(5),TFs_Multitaper_rem(6), TFs_Multitaper_rem(7),TFs_Multitaper_rem(8),TFs_Multitaper_rem(9),TFs_Multitaper_rem(10),TFs_Multitaper_rem(11),TFs_Multitaper_rem(12),TFs_Multitaper_rem(13),TFs_Multitaper_rem(14),TFs_Multitaper_rem(15), TFs_Multitaper_rem(16),TFs_Multitaper_rem(17));
GTFs_Multitaper_baseline_F = ft_freqgrandaverage(cfg,TFs_Multitaper_forg(1),TFs_Multitaper_forg(2),TFs_Multitaper_forg(3),TFs_Multitaper_forg(4),TFs_Multitaper_forg(5),TFs_Multitaper_forg(6), TFs_Multitaper_forg(7),TFs_Multitaper_forg(8),TFs_Multitaper_forg(9),TFs_Multitaper_forg(10),TFs_Multitaper_forg(11),TFs_Multitaper_forg(12),TFs_Multitaper_forg(13),TFs_Multitaper_forg(14),TFs_Multitaper_forg(15), TFs_Multitaper_forg(16),TFs_Multitaper_forg(17));
GTFs_Multitaper_baseline_K = ft_freqgrandaverage(cfg,TFs_Multitaper_k(1),TFs_Multitaper_k(2),TFs_Multitaper_k(3),TFs_Multitaper_k(4),TFs_Multitaper_k(5),TFs_Multitaper_k(6), TFs_Multitaper_k(7),TFs_Multitaper_k(8),TFs_Multitaper_k(9),TFs_Multitaper_k(10),TFs_Multitaper_k(11),TFs_Multitaper_k(12),TFs_Multitaper_k(13),TFs_Multitaper_k(14),TFs_Multitaper_k(15), TFs_Multitaper_k(16),TFs_Multitaper_k(17));
GTFs_Multitaper_baseline_KF = ft_freqgrandaverage(cfg,TFs_Multitaper_kf(1),TFs_Multitaper_kf(2),TFs_Multitaper_kf(3),TFs_Multitaper_kf(4),TFs_Multitaper_kf(5),TFs_Multitaper_kf(6), TFs_Multitaper_kf(7),TFs_Multitaper_kf(8),TFs_Multitaper_kf(9),TFs_Multitaper_kf(10),TFs_Multitaper_kf(11),TFs_Multitaper_kf(12),TFs_Multitaper_kf(13),TFs_Multitaper_kf(14),TFs_Multitaper_kf(15), TFs_Multitaper_kf(16),TFs_Multitaper_kf(17));

%% 
GTFs_diff_URUF = Ga_eR
GTFs_diff_URUF.powspctrm = Ga_eR.powspctrm - Ga_eF.powspctrm

GTFs_diff_NRNF = Ga_nR
GTFs_diff_NRNF.powspctrm = Ga_nR.powspctrm - Ga_nF.powspctrm

%% 
GTFs_diff_URUKF = Ga_eR
GTFs_diff_URUKF.powspctrm = Ga_eR.powspctrm - Ga_eKF.powspctrm

GTFs_diff_NRNKF = Ga_nR
GTFs_diff_NRNKF.powspctrm = Ga_nR.powspctrm - Ga_nKF.powspctrm

%%
GTFs_diff_URUK = Ga_eR
GTFs_diff_URUK.powspctrm = Ga_eR.powspctrm - Ga_eK.powspctrm

GTFs_diff_NRNK = Ga_eR
GTFs_diff_NRNK.powspctrm = Ga_nR.powspctrm - Ga_nK.powspctrm

%% Plot Grand Average Multitaper WITH baseline High frequency, Fig. 1d

cfg = []; 
cfg.xlim = [-.5 1.5];
cfg.ylim = [35 150];
cfg.zlim = [0 0.2];
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,Ga_eR); 
tit=title('eR amy');set(findobj(tit,'type','text'),'FontSize',12);  
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,Ga_nR);
tit=title('nR amy');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,2,3) ; ft_singleplotTFR(cfg,Ga_eKF); 
tit=title('eKF amy');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,Ga_nKF);
tit=title('nKF amy');set(findobj(tit,'type','text'),'FontSize',12);

%% Supplementary Fig 5
cfg = []; 
cfg.xlim = [-.5 1.5];
cfg.ylim = [0 34];
cfg.zlim = [0 1];
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,Ga_eR); 
tit=title('eR amy');set(findobj(tit,'type','text'),'FontSize',12);  
hold on
subplot(2,2,2); ft_singleplotTFR(cfg,Ga_nR);
tit=title('nR amy');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,2,3) ; ft_singleplotTFR(cfg,Ga_eKF); 
tit=title('eKF amy');set(findobj(tit,'type','text'),'FontSize',12);
hold on
subplot(2,2,4); ft_singleplotTFR(cfg,Ga_nKF);
tit=title('nKF amy');set(findobj(tit,'type','text'),'FontSize',12);

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

Emo_amy = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_Unpl, GTFs_Multitaper_baseline_Neu)
Mem_amy = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_R, GTFs_Multitaper_baseline_F);

EmoK_amy = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_Unplkf, GTFs_Multitaper_baseline_Neukf);
MemK_amy = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_R, GTFs_Multitaper_baseline_KF);

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

cfg.design = [ones(1,16) ones(1,16).*2;[1:16] [1:16]];

Int_amy = ft_freqstatistics(cfg,GTFs_diff_URUF, GTFs_diff_NRNF);
IntK_amy = ft_freqstatistics(cfg,GTFs_diff_URUKF, GTFs_diff_NRNKF);

%% find the range
maskt= sum(squeeze(IntK_amy.mask),1);
tvec = IntK_amy.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(IntK_amy.mask),2);
fvec = IntK_amy.freq(maskf>0);
foi=[min(fvec) max(fvec)]

%% figure 
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1.5]%'maxmin';
cfg.ylim = [35 150]% 'maxmin';
cfg.zlim = [-3 3];

E = IntK_amy;
E.powspctrm = IntK_amy.stat.*IntK_amy.mask;
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('IntK amy'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'FontSize',16) ; 
%% Fig 1e
h=figure;set(h, 'Position', get(0, 'Screensize'))

logRelative8 = IntK_amy
logRelative8.mask =IntK_amy.mask;
logRelative8.powspctrm = IntK_amy.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative8);
%%
t = toi
f = foi

pt1 = nearest(Ga_eR.time,t(1));
pt2 = nearest(Ga_eR.time,t(2));
pf1 = nearest(Ga_eR.freq,f(1));
pf2 = nearest(Ga_eR.freq,f(2));

% create matrix of conditions 
mat(:,1) = squeeze(mean(mean(Ga_eR.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(Ga_eKF.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,3) = squeeze(mean(mean(Ga_nR.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,4) = squeeze(mean(mean(Ga_nKF.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));

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
%% Fig 1f
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'));

x = [ones(16,1) ones(16,1)*2 ones(16,1)*3 ones(16,1)*4 ones(16,1)*5 ones(16,1)*6];
y = [mat(:,1) nan(16,1) mat(:,2:3) nan(16,1) mat(:,4)];
figure;beeswarm(x(:),y(:),'sort_style','up','dot_size',4,'overlay_style','sd','colormap',[1 0 0; 1 1 1;0 0 0; 0 0 1; 1 1 1;0.5 0.5 0.5])
ylim([-0.2 0.5]);
xticklabels({'eR',[],'eKF','nR',[],'nKF'})

% viy = [nan(16,1) mat(:,1)-mat(:,2) nan(16,5)];
% violinPlot(viy, 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 4,'color',  [0.9 0.9 0.9]);
% %subplot(4,7,7); hold on;
% viz = [nan(16,4) mat(:,3)-mat(:,4) nan(16,1)];
% violinPlot(viz, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 4,'color',  [0.9 0.9 0.9]);
grid off

%% STAT

cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=0;
f2=34;
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

cfg.design = [ones(1,17) ones(1,17).*2;[1:17] [1:17]];

Emo_amy_low = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_Unpl, GTFs_Multitaper_baseline_Neu)
Mem_amy_low = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_R, GTFs_Multitaper_baseline_F);
EmoK_amy_low = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_Unplkf, GTFs_Multitaper_baseline_Neukf);
MemK_amy_low = ft_freqstatistics(cfg,GTFs_Multitaper_baseline_R, GTFs_Multitaper_baseline_KF);

%% STAT

cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1.5;
f1=0;
f2=34;
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

cfg.design = [ones(1,16) ones(1,16).*2;[1:16] [1:16]];

IntK_amy_low = ft_freqstatistics(cfg,GTFs_diff_URUKF, GTFs_diff_NRNKF);


