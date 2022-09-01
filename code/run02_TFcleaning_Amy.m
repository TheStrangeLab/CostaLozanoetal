clear all;
close all;
clc;

cd C:\Users\manuela\Desktop\AversiveMemFormation\data2github\TF\Amy

%%
subjects=[2 4 6 60 13 15 16 160 25 27 32 320 21 210 10 33 34]
mlist={'s2', 's4', 's6', 's60', 's13', 's15', 's16', 's160', 's25', 's27','s32','s320','s21', 's210','sz10','s33','s34'};
behavlist = {'Patient2','Patient4','Patient6','Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient27','Patient32','Patient320','Patient21','Patient210','Patient8','PatientZ_10','Patient33','Patient34'};%

%% Do Multitaper  for Unpl Remembered

foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';%7./foi; % fRayleigh = 1/0.4s = 2.5 Hz ____ % time windows time window of 0.25 ms 0.5 ---> 1/0.5 = 2Hz
fw  = 10*ones(length(foi),1)'; %0.4*foi; %frequency smoothing
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; %1/tw(1); %4Hz tw(1) = 0.25
%cfg.foi must contain integer multiples of the Rayleigh frequency and the spectral concentration
%should also be small integer multiples of the Rayleigh frequencies, e.g.
%between 3 and 11. Cfg.tapsmofrq specifies half the spectral concentration
%spectral concentration = 2 fw(1); should be a small integer of fRayleigh  
  


for j = 1:length(subjects)
  
    load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
    
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; % length of time window
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 % time window "slides" in sec in steps of 0.01 sec (10 ms)
    cfg.pad        =    'maxperlen';  
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_E_bipolar_R(j).values = ft_freqanalysis(cfg,dataE_R);
TFs_Multitaper_HIGH_E_bipolar_R(j).values.freq = foi


end

%% Do Multitaper for Neu Remembered


foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';%7./foi; % fRayleigh = 1/0.4s = 2.5 Hz ____ % time windows time window of 0.25ms 0.5 ---> 1/0.5 = 2Hz
fw  = 10*ones(length(foi),1)'; %0.4*foi; %frequency smoothing
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; %1/tw(1); %4Hz tw(1) = 0.25
%cfg.foi must contain integer multiples of the Rayleigh frequency and the spectral concentration
%should also be small integer multiples of the Rayleigh frequencies, e.g.
%between 3 and 11. Cfg.tapsmofrq specifies half the spectral concentration
%spectral concentration = 2 fw(1); should be a small integer of fRayleigh  
  


for j = 1:length(subjects)
  
   load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
   
    %cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; % length of time window
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 % time window "slides" in sec in steps of 0.01 sec (10 ms)
    cfg.pad        =    'maxperlen'; 
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_N_bipolar_R(j).values = ft_freqanalysis(cfg,dataN_R);
TFs_Multitaper_HIGH_N_bipolar_R(j).values.freq = foi


end

%% low freq

foi = 1:1:34;

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.keeptrials = 'yes';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';

cfg.foi          = foi;
cfg.t_ftimwin    = 7./cfg.foi; 
cfg.toi          = -1.5:0.01:3;% time window "slides" in sec in steps of 0.01 sec (10 ms)
cfg.pad = 'maxperlen';

for j= 1:length(subjects)
    
        load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
        TFs_Multitaper_LOW_E_bipolar_R(j).values = ft_freqanalysis(cfg,dataE_R);
        TFs_Multitaper_LOW_E_bipolar_R(j).values.freq = foi;
      
end



for j= 1:length(subjects)
    
    load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);    
    TFs_Multitaper_LOW_N_bipolar_R(j).values = ft_freqanalysis(cfg,dataN_R);
    TFs_Multitaper_LOW_N_bipolar_R(j).values.freq = foi;

end


%% Do Multitaper  for Unpl K

foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';%7./foi; % fRayleigh = 1/0.4s = 2.5 Hz ____ % time windows time window of 0.25ms 0.5 ---> 1/0.5 = 2Hz
fw  = 10*ones(length(foi),1)'; %0.4*foi; %frequency smoothing
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; %1/tw(1); %4Hz tw(1) = 0.25
%cfg.foi must contain integer multiples of the Rayleigh frequency and the spectral concentration
%should also be small integer multiples of the Rayleigh frequencies, e.g.
%between 3 and 11. Cfg.tapsmofrq specifies half the spectral concentration
%spectral concentration = 2 fw(1); should be a small integer of fRayleigh  
  


for j = 1:length(subjects)
  
    load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
    
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; % length of time window
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 % time window "slides" in sec in steps of 0.01 sec (10 ms)
    cfg.pad        =    'maxperlen';  % hac
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_E_bipolar_K(j).values = ft_freqanalysis(cfg,dataE_K);
TFs_Multitaper_HIGH_E_bipolar_K(j).values.freq = foi


end

%% Do Multitaper for Neu Remembered


foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';%7./foi; % fRayleigh = 1/0.4s = 2.5 Hz ____ % time windows time window of 0.25ms 0.5 ---> 1/0.5 = 2Hz
fw  = 10*ones(length(foi),1)'; %0.4*foi; %frequency smoothing
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; %1/tw(1); %4Hz tw(1) = 0.25
%cfg.foi must contain integer multiples of the Rayleigh frequency and the spectral concentration
%should also be small integer multiples of the Rayleigh frequencies, e.g.
%between 3 and 11. Cfg.tapsmofrq specifies half the spectral concentration
%spectral concentration = 2 fw(1); should be a small integer of fRayleigh  
  


for j = 1:length(subjects)
  
   load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
   
    %cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; % length of time window
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 % time window "slides" in sec in steps of 0.01 sec (10 ms)
    cfg.pad        =    'maxperlen';  % hac
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_N_bipolar_K(j).values = ft_freqanalysis(cfg,dataN_K);
TFs_Multitaper_HIGH_N_bipolar_K(j).values.freq = foi


end

%% low freq

foi = 1:1:34;

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.keeptrials = 'yes';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';

cfg.foi          = foi;
cfg.t_ftimwin    = 7./cfg.foi; 
cfg.toi          = -1.5:0.01:3;% time window "slides" in sec in steps of 0.01 sec (10 ms)
cfg.pad = 'maxperlen';

for j= 1:length(subjects)
    
        load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
        TFs_Multitaper_LOW_E_bipolar_K(j).values = ft_freqanalysis(cfg,dataE_K);
        TFs_Multitaper_LOW_E_bipolar_K(j).values.freq = foi;
      
end


for j= 1:length(subjects)
    
    load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);    
    TFs_Multitaper_LOW_N_bipolar_K(j).values = ft_freqanalysis(cfg,dataN_K);
    TFs_Multitaper_LOW_N_bipolar_K(j).values.freq = foi;

end


%% Do Multitaper  for Unpl Miss

foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';%7./foi; % fRayleigh = 1/0.4s = 2.5 Hz ____ % time windows time window of 0.25ms 0.5 ---> 1/0.5 = 2Hz
fw  = 10*ones(length(foi),1)'; %0.4*foi; %frequency smoothing
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; %1/tw(1); %4Hz tw(1) = 0.25
%cfg.foi must contain integer multiples of the Rayleigh frequency and the spectral concentration
%should also be small integer multiples of the Rayleigh frequencies, e.g.
%between 3 and 11. Cfg.tapsmofrq specifies half the spectral concentration
%spectral concentration = 2 fw(1); should be a small integer of fRayleigh  
  


for j = 1:length(subjects)
  
    load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
    
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; % length of time window
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 %time window "slides" in sec in steps of 0.01 sec (10 ms)
    cfg.pad        =    'maxperlen'; 
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_E_bipolar_Miss(j).values = ft_freqanalysis(cfg,dataE_Miss);
TFs_Multitaper_HIGH_E_bipolar_Miss(j).values.freq = foi


end

%% Do Multitaper for Neu Remembered


foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';%7./foi; % fRayleigh = 1/0.4s = 2.5 Hz ____ % time windows time window of 0.25ms 0.5 ---> 1/0.5 = 2Hz
fw  = 10*ones(length(foi),1)'; %0.4*foi; %frequency smoothing
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; %1/tw(1); %4Hz tw(1) = 0.25
%cfg.foi must contain integer multiples of the Rayleigh frequency and the spectral concentration
%should also be small integer multiples of the Rayleigh frequencies, e.g.
%between 3 and 11. Cfg.tapsmofrq specifies half the spectral concentration
%spectral concentration = 2 fw(1); should be a small integer of fRayleigh  
  


for j = 1:length(subjects)
  
   load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
   
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; % length of time window
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 % time window "slides" in sec in steps of 0.01 sec (10 ms)
    cfg.pad        =    'maxperlen';  % hac
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_N_bipolar_Miss(j).values = ft_freqanalysis(cfg,dataN_Miss);
TFs_Multitaper_HIGH_N_bipolar_Miss(j).values.freq = foi


end


%% low freq

foi = 1:1:34;

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.keeptrials = 'yes';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';

cfg.foi          = foi;
cfg.t_ftimwin    = 7./cfg.foi; 
cfg.toi          = -1.5:0.01:3;% time window "slides" from -2.5 to 2.5 sec in steps of 0.01 sec (10 ms)
cfg.pad = 'maxperlen';

for j= 1:length(subjects)
    
        load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);
        TFs_Multitaper_LOW_E_bipolar_Miss(j).values = ft_freqanalysis(cfg,dataE_Miss);
        TFs_Multitaper_LOW_E_bipolar_Miss(j).values.freq = foi;
      
end



for j= 1:length(subjects)
    
    load ([behavlist{j}, '_CleanTrialsAmygdala_memory_bipolar.mat']);    
    TFs_Multitaper_LOW_N_bipolar_Miss(j).values = ft_freqanalysis(cfg,dataN_Miss);
    TFs_Multitaper_LOW_N_bipolar_Miss(j).values.freq = foi;

end

%% concatanate high low freq R

 TFs_Multitaper_E_bipolar_R = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_E_bipolar_R(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_E_bipolar_R(j).values, TFs_Multitaper_HIGH_E_bipolar_R(j).values)

 end
  
  TFs_Multitaper_N_bipolar_R = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_N_bipolar_R(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_N_bipolar_R(j).values, TFs_Multitaper_HIGH_N_bipolar_R(j).values)
 end

%% concatanate high low freq K

 TFs_Multitaper_E_bipolar_K = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_E_bipolar_K(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_E_bipolar_K(j).values, TFs_Multitaper_HIGH_E_bipolar_K(j).values)

 end
  

  TFs_Multitaper_N_bipolar_K = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_N_bipolar_K(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_N_bipolar_K(j).values, TFs_Multitaper_HIGH_N_bipolar_K(j).values)
 end


%% concatanate high low freq Miss

 TFs_Multitaper_E_bipolar_Miss = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_E_bipolar_Miss(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_E_bipolar_Miss(j).values, TFs_Multitaper_HIGH_E_bipolar_Miss(j).values)

 end
  
 
  TFs_Multitaper_N_bipolar_Miss = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_N_bipolar_Miss(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_N_bipolar_Miss(j).values, TFs_Multitaper_HIGH_N_bipolar_Miss(j).values)
 end


%% Do Baseline for emo R

TFs_Multitaper_baseline_E_bipolar_R = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_E_bipolar_R(j).values = ft_freqbaseline(cfg,TFs_Multitaper_E_bipolar_R(j).values);        
end


%% Do Baseline for neu R

TFs_Multitaper_baseline_N_bipolar_R = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_N_bipolar_R(j).values = ft_freqbaseline(cfg,TFs_Multitaper_N_bipolar_R(j).values);
end

%% Do Baseline for emo K 

TFs_Multitaper_baseline_E_bipolar_K = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_E_bipolar_K(j).values = ft_freqbaseline(cfg,TFs_Multitaper_E_bipolar_K(j).values);        
end

%% Do Baseline for neu K

TFs_Multitaper_baseline_N_bipolar_K = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_N_bipolar_K(j).values = ft_freqbaseline(cfg,TFs_Multitaper_N_bipolar_K(j).values);        
end


 %% Do Baseline for emo Miss

TFs_Multitaper_baseline_E_bipolar_Miss = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_E_bipolar_Miss(j).values = ft_freqbaseline(cfg,TFs_Multitaper_E_bipolar_Miss(j).values);        
end

%% Do Baseline for neu Miss

TFs_Multitaper_baseline_N_bipolar_Miss = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_N_bipolar_Miss(j).values = ft_freqbaseline(cfg,TFs_Multitaper_N_bipolar_Miss(j).values);
end

 %% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)% for each subject
 
    TFs_Multitaper_baseline_E_bipolar_R(j).values.rpt= [];
      
    rejectedMultitaper_E_bipolar_remembered(j).values = [];

    
    for l=1:size(TFs_Multitaper_baseline_E_bipolar_R(j).values.powspctrm,1);%each trial
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'eHit'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_eR(j) = TFs_Multitaper_baseline_E_bipolar_R(j).values;
            TFs_Multitaper_baseline_onechannel_onetrial_eR(j).powspctrm = TFs_Multitaper_baseline_E_bipolar_R(j).values.powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_eR(j).label = TFs_Multitaper_baseline_E_bipolar_R(j).values.label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eR(j))
            title(['HIGH&LOW frequencies eR, Subject:','\_', num2str(subjects(j))])


            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rejectedMultitaper_E_bipolar_remembered(j).values = [rejectedMultitaper_E_bipolar_remembered(j).values l];
            elseif artifact == 'n'
                %do nothing
            else
                return
            end
close all
        end

  
end


 %%
trials2reject = [];rpt_clean=[]; %re = [];rl = [];

for j=1:length(subjects)
    trials2reject.values = rejectedMultitaper_E_bipolar_remembered(j).values 
   
 TFs_Multitaper_baseline_E_bipolar_R(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end   

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)% for each subject
 
    TFs_Multitaper_baseline_N_bipolar_R(j).values.rpt= [];
      
    rejectedMultitaper_N_bipolar_remembered(j).values = [];

    
    for l=1:size(TFs_Multitaper_baseline_N_bipolar_R(j).values.powspctrm,1);%each trial
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'eHit'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_nR(j) = TFs_Multitaper_baseline_N_bipolar_R(j).values;
            TFs_Multitaper_baseline_onechannel_onetrial_nR(j).powspctrm = TFs_Multitaper_baseline_N_bipolar_R(j).values.powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_nR(j).label = TFs_Multitaper_baseline_N_bipolar_R(j).values.label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nR(j))
            title(['HIGH&LOW frequencies eHit, Subject:','\_', num2str(subjects(j))])


            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rejectedMultitaper_N_bipolar_remembered(j).values = [rejectedMultitaper_N_bipolar_remembered(j).values l];
            elseif artifact == 'n'
                %do nothing
            else
                return
            end
close all
        end

end


%%
trials2reject = [];rpt_clean=[]; %re = [];rl = [];

for j=1:length(subjects)
    trials2reject.values = rejectedMultitaper_N_bipolar_remembered(j).values 
   
 TFs_Multitaper_baseline_N_bipolar_R(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end   
 
 %% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)% for each subject
 
    TFs_Multitaper_baseline_E_bipolar_K(j).values.rpt= [];
      
    rejectedMultitaper_E_bipolar_familiar(j).values = [];

    
    for l=1:size(TFs_Multitaper_baseline_E_bipolar_K(j).values.powspctrm,1);%each trial
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'eK'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_eK(j) = TFs_Multitaper_baseline_E_bipolar_K(j).values;
            TFs_Multitaper_baseline_onechannel_onetrial_eK(j).powspctrm = TFs_Multitaper_baseline_E_bipolar_K(j).values.powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_eK(j).label = TFs_Multitaper_baseline_E_bipolar_K(j).values.label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eK(j))
            title(['HIGH&LOW frequencies eK, Subject:','\_', num2str(subjects(j))])


            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rejectedMultitaper_E_bipolar_familiar(j).values = [rejectedMultitaper_E_bipolar_familiar(j).values l];
            elseif artifact == 'n'
                %do nothing
            else
                return
            end
close all
        end

  
end


%%
trials2reject = [];rpt_clean=[]; %re = [];rl = [];

for j=1:length(subjects)
    trials2reject.values = rejectedMultitaper_E_bipolar_familiar(j).values 
   
 TFs_Multitaper_baseline_E_bipolar_K(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end   

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)% for each subject
 
    TFs_Multitaper_baseline_N_bipolar_K(j).values.rpt= [];
      
    rejectedMultitaper_N_bipolar_familiar(j).values = [];

    
    for l=1:size(TFs_Multitaper_baseline_N_bipolar_K(j).values.powspctrm,1);%each trial
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'nK'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_nK(j) = TFs_Multitaper_baseline_N_bipolar_K(j).values;
            TFs_Multitaper_baseline_onechannel_onetrial_nK(j).powspctrm = TFs_Multitaper_baseline_N_bipolar_K(j).values.powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_nK(j).label = TFs_Multitaper_baseline_N_bipolar_K(j).values.label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nK(j))
            title(['HIGH&LOW frequencies nK, Subject:','\_', num2str(subjects(j))])


            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rejectedMultitaper_N_bipolar_familiar(j).values = [rejectedMultitaper_N_bipolar_familiar(j).values l];
            elseif artifact == 'n'
                %do nothing
            else
                return
            end
close all
        end

  
end


%%
trials2reject = [];rpt_clean=[]; %re = [];rl = [];

for j=1:length(subjects)
    trials2reject.values = rejectedMultitaper_N_bipolar_familiar(j).values 
   
 TFs_Multitaper_baseline_N_bipolar_K(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end   

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)% for each subject
 
    TFs_Multitaper_baseline_E_bipolar_Miss(j).values.rpt= [];
      
    rejectedMultitaper_E_bipolar_forgotten(j).values = [];

    
    for l=1:size(TFs_Multitaper_baseline_E_bipolar_Miss(j).values.powspctrm,1);%each trial
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'eMiss'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j) = TFs_Multitaper_baseline_E_bipolar_Miss(j).values;
            TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j).powspctrm = TFs_Multitaper_baseline_E_bipolar_Miss(j).values.powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j).label = TFs_Multitaper_baseline_E_bipolar_Miss(j).values.label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j))
            title(['HIGH&LOW frequencies eMiss, Subject:','\_', num2str(subjects(j))])


            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rejectedMultitaper_E_bipolar_forgotten(j).values = [rejectedMultitaper_E_bipolar_forgotten(j).values l];
           
            elseif artifact == 'n'
                %do nothing
            else
                return
            end
close all
        end

  
end


%%
trials2reject = [];rpt_clean=[]; %re = [];rl = [];

for j=1:length(subjects)
    trials2reject.values = rejectedMultitaper_E_bipolar_forgotten(j).values 
   
 TFs_Multitaper_baseline_E_bipolar_Miss(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end   

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)% for each subject
 
    TFs_Multitaper_baseline_N_bipolar_Miss(j).values.rpt= [];
      
    rejectedMultitaper_N_bipolar_forgotten(j).values = [];

    
    for l=1:size(TFs_Multitaper_baseline_N_bipolar_Miss(j).values.powspctrm,1);%each trial
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'nMiss'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j) = TFs_Multitaper_baseline_N_bipolar_Miss(j).values;
            TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j).powspctrm = TFs_Multitaper_baseline_N_bipolar_Miss(j).values.powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j).label = TFs_Multitaper_baseline_N_bipolar_Miss(j).values.label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j))
            title(['HIGH&LOW frequencies nMiss, Subject:','\_', num2str(subjects(j))])


            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rejectedMultitaper_N_bipolar_forgotten(j).values = [rejectedMultitaper_N_bipolar_forgotten(j).values l];
             
            elseif artifact == 'n'
                %do nothing
            else
                return
            end
close all
        end

  
end

%%
trials2reject = [];rpt_clean=[]; %re = [];rl = [];

for j=1:length(subjects)
    trials2reject.values = rejectedMultitaper_N_bipolar_forgotten(j).values 
   
 TFs_Multitaper_baseline_N_bipolar_Miss(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end   

%% 
eR = TFs_Multitaper_baseline_eR
nR = TFs_Multitaper_baseline_nR
eK = TFs_Multitaper_baseline_eK
nK = TFs_Multitaper_baseline_nK
eF = TFs_Multitaper_baseline_eMiss
nF = TFs_Multitaper_baseline_nMiss
%%
for j = 1:length(subjects)
    cfg= []
    cfg.parameter = 'powspctrm'
    eKF(j).values = ft_appendfreq(cfg, eK(j).values, eF(j).values);
    nKF(j).values = ft_appendfreq(cfg, nK(j).values, nF(j).values);
end

% save(TimeFreq_cond_allamysubjects.mat, 'eR','eK','eF','nR','nK','nF','eKF','nKF');