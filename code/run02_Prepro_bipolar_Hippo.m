clear all
close all
clc


%% make subject list and define channels to use

subjects= [60 13 15 16 160 25 32 33]
mlist={'s60', 's13', 's15', 's16', 's160', 's25', 's32', 's33'} 
behavlist = {'Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient32','Patient33'}
        
%%s60r
chanofinterest{1} = {'THD1', 'THD2'};
labelorg{1} = {'THD1', 'THD2'};
bipolabel{1} = {'THD1 - THD2'};
bipolmat{1} =  [+1 -1];  
%%s13
chanofinterest{2} = {'HA2', 'HA3', 'HA4'};
labelorg{2} = {'HA2', 'HA3', 'HA4'};
bipolabel{2} = {'HA2-HA3','HA3-HA4'};
bipolmat{2} = [ +1 -1 0; 0 +1 -1];   
%%s15
chanofinterest{3} = {'THA2', 'THA3'};
labelorg{3} = {'THA2', 'THA3'};
bipolabel{3} = {'THA2 - THA3'};
bipolmat{3} = [+1 -1];   
%%s16
chanofinterest{4} = {'HI2', 'HI3', 'HI4'};
labelorg{4} = {'HI2', 'HI3', 'HI4'};
bipolabel{4} = {'HI2 - HI3', 'HI3 - HI4'};
bipolmat{4} = [ +1 -1 0; 0 +1 -1];      
%%s16r 
chanofinterest{5} = {'HD2', 'HD3', 'HD4'};
labelorg{5} = {'HD2', 'HD3', 'HD4'};
bipolabel{5} = {'HD2 - HD3','HD3 - HD4'};
bipolmat{5} = [ +1 -1 0; 0 +1 -1];    
%s25
chanofinterest{6} = {'HI3', 'HI4'};
labelorg{6} = {'HI3', 'HI4'};
bipolabel{6} = {'HI3 - HI4'};
bipolmat{6} = [ +1 -1];    
%s32
chanofinterest{7} = {'D1', 'D2'};
labelorg{7} = {'D1', 'D2'};
bipolabel{7} = {'D1 - D2'};
bipolmat{7} = [+1 -1]; 
%s33
chanofinterest{8} = {'TB1', 'TB2'};
labelorg{8} = {'TB1', 'TB2'};
bipolabel{8} = {'TB1 - TB2'};
bipolmat{8} = [+1 -1];   
%%
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

    
    % get eR trials only
    trl = cfg.trl; %remember all trials
    cfg.trl = trl(enc_onsets.eCorrRem,:);
    

    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};

    dataeR_bipolar = ft_preprocessing(cfg);
  
    
   % get nR trials only
   
    cfg.trl = trl(enc_onsets.nCorrRem,:);
    
    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};

    datanR_bipolar = ft_preprocessing(cfg);
    
    % get eK trials only
    cfg.trl = trl(enc_onsets.eCorrFam,:);

    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};

    dataeK_bipolar = ft_preprocessing(cfg);
  
    
   % get nK trials only
   cfg.trl = trl(enc_onsets.nCorrFam,:);
 
    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};

    datanK_bipolar = ft_preprocessing(cfg);  
    
    % get eMiss trials only
    cfg.trl = trl(enc_onsets.eMissed,:);
   
    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};

    dataeMiss_bipolar = ft_preprocessing(cfg);
  
    
   % get nK trials only
    cfg.trl = trl(enc_onsets.nMissed,:);

    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};

    datanMiss_bipolar = ft_preprocessing(cfg);  
    
% look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eRclean = ft_rejectvisual(cfg,dataeR_bipolar);
    
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nRclean = ft_rejectvisual(cfg,datanR_bipolar);
    
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eFclean = ft_rejectvisual(cfg,dataeMiss_bipolar);
    
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nFclean = ft_rejectvisual(cfg,datanMiss_bipolar);

    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eKclean = ft_rejectvisual(cfg,dataeK_bipolar);
     
   % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nKclean = ft_rejectvisual(cfg,datanK_bipolar);

    %% save preprocessed data

    %file = [behavlist{v},'_CleanTrialsHippo_Memory_bipolarL.mat'];
    %save(file,'dataeR_bipolar','datanR_bipolar','dataeK_bipolar','datanK_bipolar','dataeMiss_bipolar','datanMiss_bipolar');

end