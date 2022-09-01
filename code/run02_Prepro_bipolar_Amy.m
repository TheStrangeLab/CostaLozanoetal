clear all
close all
clc

%% make subject list and define channels to use

subjects=[2 4 6 60 13 15 16 160 25 27 32 320 21 210 33 34]
mlist={'s2', 's4', 's6', 's60', 's13', 's15', 's16', 's160', 's25', 's27','s32','s320','s21', 's210','s33','s34'};
behavlist = {'Patient2','Patient4','Patient6','Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient27','Patient32','Patient320','Patient21','Patient210','Patient8','Patient33','Patient34'};%

%s2
chanofinterest{1} = {'TA1','TA2','TA3'};
labelorg{1} = {'TA1','TA2','TA3'};
bipolabel{1} = {'TA1 - TA2','TA2 - TA3'};
bipolmat{1} = [+1 -1 0; 0 +1 -1];
%%s4
chanofinterest{2} = {'TA1','TA2','TA3','TA4'};
labelorg{2} = {'TA1','TA2','TA3','TA4'};
bipolabel{2} = {'TA1 - TA2','TA2 - TA3','TA3 -TA4'};
bipolmat{2} = [ +1 -1 0 0; 0 +1 -1 0; 0 0 +1 -1];
%%s6
chanofinterest{3} = {'TAI1','TAI2','TAI3'};
labelorg{3} = {'TAI1','TAI2','TAI3'};
bipolabel{3} = {'TAI1-TAI2','TAI2-TAI3'};
bipolmat{3} = [ +1 -1 0; 0 +1 -1];
%%s6r
chanofinterest{4} = {'TAD1','TAD2','TAD3'};
labelorg{4} = {'TAD1','TAD2','TAD3'};
bipolabel{4} = {'TAD1-TAD2','TAD2-TAD3'};
bipolmat{4} =  [ +1 -1 0; 0 +1 -1];
%%s13
chanofinterest{5} = {'AM1','AM2','AM3'};
labelorg{5} = {'AM1','AM2','AM3'};
bipolabel{5} = {'AM1 -AM2','AM2 -AM3'};
bipolmat{5} = [ +1 -1 0; 0 +1 -1];
%%s15
chanofinterest{6} = {'A1','A2','A3'};
labelorg{6} = {'A1','A2','A3'};
bipolabel{6} = {'A1-A2','A2-A3'};
bipolmat{6} = [ +1 -1 0; 0 +1 -1];
%%s16
chanofinterest{7} = {'AI1','AI2'};
labelorg{7} = {'AI1','AI2'};
bipolabel{7} = {'AI1 - AI2'};
bipolmat{7} = [ +1 -1];
%%s16r
chanofinterest{8} = {'AD1','AD2'};
labelorg{8} = {'AD1','AD2'};
bipolabel{8} = {'AD1 - AD2'};
bipolmat{8} = [ +1 -1];
%%s25
chanofinterest{9} = {'A3','A4'};
labelorg{9} = {'A3','A4'};
bipolabel{9} = {'A3 -A4' };
bipolmat{9} = [ +1 -1];
%%s27
chanofinterest{10} = {'T1A1','T1A2'};
labelorg{10} = {'T1A1','T1A2'};
bipolabel{10} = {'T1A1-T1A2'};
bipolmat{10} = [-1 +1];
%%s32
chanofinterest{11} = {'B2','B3'};
labelorg{11} = {'B2','B3'};
bipolabel{11} = {'B2 - B3'};
bipolmat{11} = [ +1 -1];
%%s320
chanofinterest{12} = {'F2','F3'};
labelorg{12} = {'F2','F3'};
bipolabel{12} = {'F2 - F3'};
bipolmat{12} = [ +1 -1];
%%s21
chanofinterest{13} = {'AI1','AI2','AI3'};
labelorg{13} = {'AI1','AI2','AI3'};
bipolabel{13} = {'AI1 - AI2','AI2 - AI3'};
bipolmat{13} = [ +1 -1 0; 0 +1 -1];
%%s210
chanofinterest{14} = {'AD1','AD2',};
labelorg{14} = {'AD1','AD2',};
bipolabel{14} = {'AD1 - AD2'};
bipolmat{14} = [ +1 -1];
%%s33
chanofinterest{15} = {'A1', 'A2', 'A3'};
labelorg{15} = {'A1', 'A2', 'A3'};
bipolabel{15} = {'A1 - A2', 'A2 - A3'};
bipolmat{15} = [+1 -1 0; 0 +1 -1];
%%s34
chanofinterest{16} = {'B2', 'B3', 'B4'};
labelorg{16} = {'B2', 'B3', 'B4'};
bipolabel{16} = {'B2 - B3', 'B3 - B4'};
bipolmat{16} = [+1 -1 0; 0 +1 -1];
%%

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))

v=0;
for sub=subjects
    %% define subject
    v=v+1;
    sdir=['~\DataRuber\Patient',num2str(sub)];
    cd(sdir)
    
    cd Enc
    cond = load('IAPS_enc.mat')
    
    cd(sdir)
    
    filename = [mlist{v},'_IAPS_enc.edf'];
    
    % define trials
    cfg=[];
    cfg.dataset= filename;
    cfg.prestim = 7.5;
    cfg.poststim = 7.5;
    cfg.trialfun = 'mytrialfun_IntraCranial_iaps';
    cfg.triggers = cond.IAPS_enc;
    cfg = ft_definetrial(cfg);
    
    % load behavioral data
    load('enc_onsets.mat');
    
    
    % get eR trials only
    trl = cfg.trl; %remember all trials
    cfg.trl = trl(enc_onsets.eCorrRem,:);
    
    % now read the data
    
    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};
    
    dataE_R = ft_preprocessing(cfg);
    
    
    % get nR trials only
    
    cfg.trl = trl(enc_onsets.nCorrRem,:);
    
    % now read the data
    
    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};
    
    dataN_R = ft_preprocessing(cfg);
    
    % get eK trials only
    cfg.trl = trl(enc_onsets.eCorrFam,:);
    
    
    % now read the data
    
    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};
    
    dataE_K = ft_preprocessing(cfg);
    
    
    % get nK trials only
    cfg.trl = trl(enc_onsets.nCorrFam,:);
    
    % now read the data
    
    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};
    
    dataN_K = ft_preprocessing(cfg);
    
    % get eMiss trials only
    cfg.trl = trl(enc_onsets.eMissed,:);
    
    
    % now read the data
    
    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};
    
    dataE_Miss = ft_preprocessing(cfg);
    
    
    % get nK trials only
    cfg.trl = trl(enc_onsets.nMissed,:);
    
    % now read the data
    
    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};
    
    dataN_Miss = ft_preprocessing(cfg);
    
    % look for spikes
        cfg          = [];
        cfg.channel = bipolabel{v};
        cfg.latency = [-0.5 1.5]
        cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
        cfg.method   = 'channel';
        data_eRclean = ft_rejectvisual(cfg,dataE_R);
    
        % look for spikes
        cfg          = [];
        cfg.channel = bipolabel{v};
        cfg.latency = [-0.5 1.5]
        cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
        cfg.method   = 'channel';
        data_nRclean = ft_rejectvisual(cfg,dataN_R);
    
        % look for spikes
        cfg          = [];
        cfg.channel = bipolabel{v};
        cfg.latency = [-0.5 1.5]
        cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
        cfg.method   = 'channel';
        data_eFclean = ft_rejectvisual(cfg,dataE_Miss);
    
        % look for spikes
        cfg          = [];
        cfg.channel = bipolabel{v};
        cfg.latency = [-0.5 1.5]
        cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
        cfg.method   = 'channel';
        data_nFclean = ft_rejectvisual(cfg,dataN_Miss);
    
        % look for spikes
        cfg          = [];
        cfg.channel = bipolabel{v};
        cfg.latency = [-0.5 1.5]
        cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
        cfg.method   = 'channel';
        data_eKclean = ft_rejectvisual(cfg,dataE_K);
    
       % look for spikes
        cfg          = [];
        cfg.channel = bipolabel{v};
        cfg.latency = [-0.5 1.5]
        cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
        cfg.method   = 'channel';
        data_nKclean = ft_rejectvisual(cfg,dataN_K);
    
    %% save preprocessed data
    
    %file = [behavlist{v},'_CleanTrialsAmygdala_memory_bipolar.mat'];
    %save(file,'dataE_R','dataN_R','dataE_K','dataN_K','dataE_Miss','dataN_Miss');
       
end