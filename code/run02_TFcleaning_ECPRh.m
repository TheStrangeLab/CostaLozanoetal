clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

cd C:\Users\manuela\Desktop\AversiveMemFormation\data2github\TF\ECPRh

%%
subjects=[13 15 34 5 500 600 8 10];
mlist={'s13', 's15','s34','s5','s500','s600','s8','s10'};
behavlist = {'Patient13','Patient15','Patient34','Patientz5','Patientz500','Patientz600','Patientz8','Patientz10'};
%%
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
  
j=0;
for j=1:length(subjects)
    
    load ([behavlist{j}, '_dws_ParaHcECPC_MemoryK.mat']);
 
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; % length of time window
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 %-2.5:0.01:2.5; % time window "slides" in sec in steps of 0.01 sec (10 ms)
    cfg.pad        =    'maxperlen';  % hac
    cfg.taper   = 'dpss'

TFs_Multitaper_HIGH_eR(j) = ft_freqanalysis(cfg,data_eRclean);
TFs_Multitaper_HIGH_eR(j).freq = foi

TFs_Multitaper_HIGH_nR(j) = ft_freqanalysis(cfg,data_nRclean);
TFs_Multitaper_HIGH_nR(j).freq = foi

TFs_Multitaper_HIGH_eF(j) = ft_freqanalysis(cfg,data_eFclean);
TFs_Multitaper_HIGH_eF(j).freq = foi

TFs_Multitaper_HIGH_nF(j) = ft_freqanalysis(cfg,data_nFclean);
TFs_Multitaper_HIGH_nF(j).freq = foi

TFs_Multitaper_HIGH_eK(j) = ft_freqanalysis(cfg,data_eKclean);
TFs_Multitaper_HIGH_eK(j).freq = foi

TFs_Multitaper_HIGH_nK(j) = ft_freqanalysis(cfg,data_nKclean);
TFs_Multitaper_HIGH_nK(j).freq = foi

end


%%
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

j=0;
for j=1:length(subjects)
    
    load ([behavlist{j}, '_dws_ParaHcECPC_MemoryK.mat']);
    
    TFs_Multitaper_LOW_eR(j) = ft_freqanalysis(cfg,data_eRclean);
    TFs_Multitaper_LOW_eR(j).freq = foi;
        
    TFs_Multitaper_LOW_nR(j) = ft_freqanalysis(cfg,data_nRclean);
    TFs_Multitaper_LOW_nR(j).freq = foi;
    
    TFs_Multitaper_LOW_eF(j) = ft_freqanalysis(cfg,data_eFclean);
    TFs_Multitaper_LOW_eF(j).freq = foi;
        
    TFs_Multitaper_LOW_nF(j) = ft_freqanalysis(cfg,data_nFclean);
    TFs_Multitaper_LOW_nF(j).freq = foi;
    
    TFs_Multitaper_LOW_eK(j) = ft_freqanalysis(cfg,data_eKclean);
    TFs_Multitaper_LOW_eK(j).freq = foi;
    
    TFs_Multitaper_LOW_nK(j) = ft_freqanalysis(cfg,data_nKclean);
    TFs_Multitaper_LOW_nK(j).freq = foi;
      
end
%% concatanate high low freq
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_eR(j) = ft_appendfreq(cfg, TFs_Multitaper_LOW_eR(j), TFs_Multitaper_HIGH_eR(j))
        TFs_Multitaper_nR(j) = ft_appendfreq(cfg, TFs_Multitaper_LOW_nR(j), TFs_Multitaper_HIGH_nR(j))
        TFs_Multitaper_eF(j) = ft_appendfreq(cfg, TFs_Multitaper_LOW_eF(j), TFs_Multitaper_HIGH_eF(j))
        TFs_Multitaper_nF(j) = ft_appendfreq(cfg, TFs_Multitaper_LOW_nF(j), TFs_Multitaper_HIGH_nF(j))
        TFs_Multitaper_eK(j) = ft_appendfreq(cfg, TFs_Multitaper_LOW_eK(j), TFs_Multitaper_HIGH_eK(j))
        TFs_Multitaper_nK(j) = ft_appendfreq(cfg, TFs_Multitaper_LOW_nK(j), TFs_Multitaper_HIGH_nK(j))
 end

%% Do Baseline

for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_eR(j) = ft_freqbaseline(cfg,TFs_Multitaper_eR(j));
        TFs_Multitaper_baseline_nR(j) = ft_freqbaseline(cfg,TFs_Multitaper_nR(j));
        TFs_Multitaper_baseline_eF(j) = ft_freqbaseline(cfg,TFs_Multitaper_eF(j));
        TFs_Multitaper_baseline_nF(j) = ft_freqbaseline(cfg,TFs_Multitaper_nF(j));
        TFs_Multitaper_baseline_eK(j) = ft_freqbaseline(cfg,TFs_Multitaper_eK(j));
        TFs_Multitaper_baseline_nK(j) = ft_freqbaseline(cfg,TFs_Multitaper_nK(j));
end

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 1.5];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_eR(j).rpt= [];
      
    rj_er(j).values = [];
    
    for l=1:size(TFs_Multitaper_baseline_eR(j).powspctrm,1);
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'eR'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_eR(j) = TFs_Multitaper_baseline_eR(j);
            TFs_Multitaper_baseline_onechannel_onetrial_eR(j).powspctrm = TFs_Multitaper_baseline_eR(j).powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_eR(j).label = TFs_Multitaper_baseline_eR(j).label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eR(j))
            title(['HIGH&LOW frequencies eR, Subject:','\_', num2str(subjects(j))])
            
            pause(1) % you can change this depending on how long do you want to see the plot
            close all 
            
            
            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rj_er(j).values = [rj_er(j).values l];
            elseif artifact == 'n'
            else
                return
            end
    end
end

%%
trials2reject = [];rpt_clean=[];

for j=1:length(subjects)
trials2reject = rj_er(j).values 
   
TFs_Multitaper_baseline_eR(j).powspctrm(trials2reject, :,:,:) = [];

end   

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 1.5];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_eK(j).rpt= [];
      
    rj_ek(j).values = [];
    
    for l=1:size(TFs_Multitaper_baseline_eK(j).powspctrm,1);
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'eK'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_eK(j) = TFs_Multitaper_baseline_eK(j);
            TFs_Multitaper_baseline_onechannel_onetrial_eK(j).powspctrm = TFs_Multitaper_baseline_eK(j).powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_eK(j).label = TFs_Multitaper_baseline_eK(j).label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eK(j))
            title(['HIGH&LOW frequencies eK, Subject:','\_', num2str(subjects(j))])
            
            pause(1) 
            close all 
            
            
            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rj_ek(j).values = [rj_ek(j).values l];
            elseif artifact == 'n'
            else
                return
            end
    end
end

%%
trials2reject = [];rpt_clean=[];

for j=1:length(subjects)
trials2reject = rj_ek(j).values 
   
TFs_Multitaper_baseline_eK(j).powspctrm(trials2reject, :,:,:) = [];

end  

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 1.5];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_eF(j).rpt= [];
      
    rj_ef(j).values = [];
    
    for l=1:size(TFs_Multitaper_baseline_eF(j).powspctrm,1);
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'eF'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_eF(j) = TFs_Multitaper_baseline_eF(j);
            TFs_Multitaper_baseline_onechannel_onetrial_eF(j).powspctrm = TFs_Multitaper_baseline_eF(j).powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_eF(j).label = TFs_Multitaper_baseline_eF(j).label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eF(j))
            title(['HIGH&LOW frequencies eF, Subject:','\_', num2str(subjects(j))])
            
            pause(1) 
            close all 
            
            
            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rj_ef(j).values = [rj_ef(j).values l];
            elseif artifact == 'n'
            else
                return
            end
    end
end

%%
trials2reject = [];rpt_clean=[];

for j=1:length(subjects)
trials2reject = rj_ef(j).values 
   
TFs_Multitaper_baseline_eF(j).powspctrm(trials2reject, :,:,:) = [];
end  

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 1.5];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_nR(j).rpt= [];
      
    rj_nr(j).values = [];
    
    for l=1:size(TFs_Multitaper_baseline_nR(j).powspctrm,1);
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'nR'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_nR(j) = TFs_Multitaper_baseline_nR(j);
            TFs_Multitaper_baseline_onechannel_onetrial_nR(j).powspctrm = TFs_Multitaper_baseline_nR(j).powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_nR(j).label = TFs_Multitaper_baseline_nR(j).label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nR(j))
            title(['HIGH&LOW frequencies nR, Subject:','\_', num2str(subjects(j))])
            
            close all 
            
            
            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rj_nr(j).values = [rj_nr(j).values l];
            elseif artifact == 'n'
            else
                return
            end
    end
end

%%
trials2reject = [];rpt_clean=[];

for j=1:length(subjects)
trials2reject = rj_nr(j).values 
   
TFs_Multitaper_baseline_nR(j).powspctrm(trials2reject, :,:,:) = [];

end   

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 1.5];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_nK(j).rpt= [];
      
    rj_nk(j).values = [];
    
    for l=1:size(TFs_Multitaper_baseline_nK(j).powspctrm,1);
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'nK'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_nK(j) = TFs_Multitaper_baseline_nK(j);
            TFs_Multitaper_baseline_onechannel_onetrial_nK(j).powspctrm = TFs_Multitaper_baseline_nK(j).powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_nK(j).label = TFs_Multitaper_baseline_nK(j).label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nK(j))
            title(['HIGH&LOW frequencies nK, Subject:','\_', num2str(subjects(j))])
            
            close all 
            
            
            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rj_nk(j).values = [rj_nk(j).values l];
            elseif artifact == 'n'
            else
                return
            end
    end
end

%%
trials2reject = [];rpt_clean=[];

for j=1:length(subjects)
trials2reject = rj_nk(j).values 
   
TFs_Multitaper_baseline_nK(j).powspctrm(trials2reject, :,:,:) = [];

end  

%% Artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 1.5];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_nF(j).rpt= [];
      
    rj_nf(j).values = [];
    
    for l=1:size(TFs_Multitaper_baseline_nF(j).powspctrm,1);
           
            h = figure;
            set (h, 'Units', 'normalized', 'Position', [0,0,1,1]);
            switch i
                case 1
                    label = 'nF'
            end

            %take only the channel you are interested in
            TFs_Multitaper_baseline_onechannel_onetrial_nF(j) = TFs_Multitaper_baseline_nF(j);
            TFs_Multitaper_baseline_onechannel_onetrial_nF(j).powspctrm = TFs_Multitaper_baseline_nF(j).powspctrm(l,:,:,:);
            TFs_Multitaper_baseline_onechannel_onetrial_nF(j).label = TFs_Multitaper_baseline_nF(j).label;

%             subplot(2,1,1);
            ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nF(j))
            title(['HIGH&LOW frequencies nF, Subject:','\_', num2str(subjects(j))])
            
            pause(1) 
            close all 
            
            
            fprintf(strcat('******************************************trial:_',num2str(l)))
            fprintf('\n')
            artifact = 'a'; %initialization different to y or n to enter the loop
            while (artifact ~= 'y' && artifact ~= 'n')
                artifact = input('Is this trial an artifact? (y/n): ','s');
            end
            if artifact == 'y'
                rj_nf(j).values = [rj_nf(j).values l];
            elseif artifact == 'n'
            else
                return
            end
    end
end

%%
trials2reject = [];rpt_clean=[];

for j=1:length(subjects)
trials2reject = rj_nf(j).values 
   
TFs_Multitaper_baseline_nF(j).powspctrm(trials2reject, :,:,:) = [];

end 

%% Do Avaragae with baseline 
%collapsing channels All
TFs_Multitaper_baseline_average_eR = [];
for j = 1:length(subjects)
   
        TFs_Multitaper_baseline_average_eR(j).values  = TFs_Multitaper_baseline_eR(j);
        if size(TFs_Multitaper_baseline_eR(j).powspctrm,2) ~= 1 %channels
            TFs_Multitaper_baseline_average_eR(j).values .powspctrm = mean(TFs_Multitaper_baseline_eR(j).powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TFs_Multitaper_baseline_average_eR(j).values .label = {'ParaHcEPC'};
        size(TFs_Multitaper_baseline_average_eR(j).values.powspctrm)
    end

%% Do Average WITH Baseline correction for Neutral R

%collapsing channels All
TFs_Multitaper_baseline_average_nR = [];
for j = 1:length(subjects)
   
        TFs_Multitaper_baseline_average_nR(j).values = TFs_Multitaper_baseline_nR(j);
        if size(TFs_Multitaper_baseline_nR(j).powspctrm,2) ~= 1 %channels
            TFs_Multitaper_baseline_average_nR(j).values.powspctrm = mean(TFs_Multitaper_baseline_nR(j).powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TFs_Multitaper_baseline_average_nR(j).values.label = {'ParaHcEPC'};
        size(TFs_Multitaper_baseline_average_nR(j).values.powspctrm)
    end

%%
TFs_Multitaper_baseline_average_eF = [];
for j = 1:length(subjects)
   
        TFs_Multitaper_baseline_average_eF(j).values = TFs_Multitaper_baseline_eF(j);
        if size(TFs_Multitaper_baseline_eF(j).powspctrm,2) ~= 1 %channels
            TFs_Multitaper_baseline_average_eF(j).values.powspctrm = mean(TFs_Multitaper_baseline_eF(j).powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TFs_Multitaper_baseline_average_eF(j).values.label = {'ParaHcEPC'};
        size(TFs_Multitaper_baseline_average_eF(j).values.powspctrm)
end

%%
TFs_Multitaper_baseline_average_eK = [];
for j = 1:length(subjects)
   
        TFs_Multitaper_baseline_average_eK(j).values = TFs_Multitaper_baseline_eK(j);
        if size(TFs_Multitaper_baseline_eK(j).powspctrm,2) ~= 1 %channels
            TFs_Multitaper_baseline_average_eK(j).values.powspctrm = mean(TFs_Multitaper_baseline_eK(j).powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TFs_Multitaper_baseline_average_eK(j).values.label = {'ParaHcEPC'};
        size(TFs_Multitaper_baseline_average_eK(j).values.powspctrm)
end

%% Do Average WITH Baseline correction for Neutral R

%collapsing channels All
TFs_Multitaper_baseline_average_nF = [];
for j = 1:length(subjects)
   
        TFs_Multitaper_baseline_average_nF(j).values = TFs_Multitaper_baseline_nF(j);
        if size(TFs_Multitaper_baseline_nF(j).powspctrm,2) ~= 1 %channels
            TFs_Multitaper_baseline_average_nF(j).values.powspctrm = mean(TFs_Multitaper_baseline_nF(j).powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TFs_Multitaper_baseline_average_nF(j).values.label = {'ParaHcEPC'};
        size(TFs_Multitaper_baseline_average_nF(j).values.powspctrm)
end

%%
TFs_Multitaper_baseline_average_nK = [];
for j = 1:length(subjects)
   
        TFs_Multitaper_baseline_average_nK(j).values = TFs_Multitaper_baseline_nK(j);
        if size(TFs_Multitaper_baseline_nK(j).powspctrm,2) ~= 1 %channels
            TFs_Multitaper_baseline_average_nK(j).values.powspctrm = mean(TFs_Multitaper_baseline_nK(j).powspctrm,2); % I want to mean the channel큦 dimension 
        end
        TFs_Multitaper_baseline_average_nK(j).values.label = {'ParaHcEPC'};
        size(TFs_Multitaper_baseline_average_nK(j).values.powspctrm)
end

% cd C:\Users\manuela\Desktop\AversiveMemFormation\data2github\TF\ECPRh
% save(TimeFreq_cond_ECPRhsubjects.mat, 'TFs_Multitaper_baseline_average_eR','TFs_Multitaper_baseline_average_eK','TFs_Multitaper_baseline_average_eF','TFs_Multitaper_baseline_average_nR','TFs_Multitaper_baseline_average_nK','TFs_Multitaper_baseline_average_nF');