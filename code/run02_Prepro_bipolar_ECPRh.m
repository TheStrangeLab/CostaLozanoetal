clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

%% make subject list

subjects=[13 15 34];
mlist={'s13', 's15', 's34'};
behavlist = {'Patient13','Patient15','Patient34'};
   
%%s13 - entorhinal
chanofinterest{1} = {'TBA1','TBA2'};
labelorg{1} = {'TBA1','TBA2'};
bipolabel{1} = {'TBA1 - TBA2'};
bipolmat{1} = [+1 -1];             

%%s15 - perirhinal/entorhinal
chanofinterest{2} = {'TBA1','TBA2'};
labelorg{2} = {'TBA1','TBA2'};
bipolabel{2} = {'TBA1 - TBA2'};
bipolmat{2} = [+1 -1]; 

%%s34 -perirhinal/entorhinal
chanofinterest{3} = {'D1','D2'};
labelorg{3} = {'D1','D2'};
bipolabel{3} = {'D1 - D2'};
bipolmat{3} = [+1 -1];
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
    cfg.prestim = 7.5; %1.5;
    cfg.poststim = 7.5; %3; 
    cfg.trialfun = 'mytrialfun_IntraCranial_iaps';
    cfg.triggers = cond.IAPS_enc;
    cfg = ft_definetrial(cfg);

    
    % load behavioral data
    load('enc_onsets.mat');
   
            
    % get emotional trials only
    trl = cfg.trl; %remember all trials     
    cfg.trl = trl(enc_onsets.eCorrRem,:);

    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};

    data_eR = ft_preprocessing(cfg);
  
    cfg.trl = trl(enc_onsets.eMissed,:);
    data_eF = ft_preprocessing(cfg);
    
    
    cfg.trl = trl(enc_onsets.nCorrRem,:);
    data_nR = ft_preprocessing(cfg);
    
    
    cfg.trl = trl(enc_onsets.nCorrFam,:);
    data_nK = ft_preprocessing(cfg);
    
    cfg.trl = trl(enc_onsets.nMissed,:);
    data_nF = ft_preprocessing(cfg);
    
    cfg.trl = trl(enc_onsets.eCorrFam,:);
    data_eK = ft_preprocessing(cfg);
    
    
   % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eRclean = ft_rejectvisual(cfg,data_eR);
    
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nRclean = ft_rejectvisual(cfg,data_nR);
    
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eFclean = ft_rejectvisual(cfg,data_eF);
    
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nFclean = ft_rejectvisual(cfg,data_nF);

    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eKclean = ft_rejectvisual(cfg,data_eK);
     
   % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nKclean = ft_rejectvisual(cfg,data_nK);
  %% save preprocessed data   
  

    file = [behavlist{v},'_CleanTrialsParaHcECPC_memory_bipolar.mat'];
save(file,'data_eRclean','data_eKclean','data_eFclean','data_nRclean','data_nKclean','data_nFclean');
  

end
