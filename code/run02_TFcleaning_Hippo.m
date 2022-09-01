clear all;
close all;
clc;

cd C:\Users\manuela\Desktop\AversiveMemFormation\data2github\TF\Hippo
%%
subjects=[60 13 15 16 160 25 32 10 33];
mlist={'s60', 's13', 's15', 's16', 's160', 's25','s32','s10','s33'}
behavlist = {'Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient32','PatientZ_10', 'Patient33'};%

%% Do Multitaper for R

foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';
fw  = 10*ones(length(foi),1)'; 
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw;
  
for j = 1:length(subjects)
  
    load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
    
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; 
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 
    cfg.pad        =    'maxperlen'; 
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_eR(j).values = ft_freqanalysis(cfg,dataeR_bipolar);
TFs_Multitaper_HIGH_eR(j).values.freq = foi
end


foi = 35:2.5:150; % Frequency of interest
tw  = 0.4*ones(length(foi),1)';
fw  = 10*ones(length(foi),1)'; 
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw;   
  


for j = 1:length(subjects)
  
   load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
   
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; 
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3 
    cfg.pad        =    'maxperlen'; 
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_nR(j).values = ft_freqanalysis(cfg,datanR_bipolar);
TFs_Multitaper_HIGH_nR(j).values.freq = foi


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
cfg.toi          = -1.5:0.01:3;
cfg.pad = 'maxperlen';
   
for j= 1:length(subjects)
    
        load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
        TFs_Multitaper_LOW_eR(j).values = ft_freqanalysis(cfg,dataeR_bipolar);
        TFs_Multitaper_LOW_eR(j).values.freq = foi;
      
end



for j= 1:length(subjects)
    
    load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);    
    TFs_Multitaper_LOW_nR(j).values = ft_freqanalysis(cfg,datanR_bipolar);
    TFs_Multitaper_LOW_nR(j).values.freq = foi;

end


%% Do Multitaper for K

foi = 35:2.5:150; 
tw  = 0.4*ones(length(foi),1)';
fw  = 10*ones(length(foi),1)';
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; 
  


for j = 1:length(subjects)
   load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; 
    cfg.tapsmofrq = fw;
    cfg.toi        = -1.5:0.01:3
    cfg.pad        =    'maxperlen'; 
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_eK(j).values = ft_freqanalysis(cfg,dataeK_bipolar);
TFs_Multitaper_HIGH_eK(j).values.freq = foi


end


for j= 1:length(subjects)
     load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
    TFs_Multitaper_HIGH_nK(j).values = ft_freqanalysis(cfg,datanK_bipolar);
    TFs_Multitaper_HIGH_nK(j).values.freq = foi;

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
cfg.toi          = -1.5:0.01:3;
cfg.pad = 'maxperlen';
   
for j= 1:length(subjects)
    load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
    
        TFs_Multitaper_LOW_eK(j).values = ft_freqanalysis(cfg,dataeK_bipolar);
        TFs_Multitaper_LOW_eK(j).values.freq = foi;
      
end


for j= 1:length(subjects)
    load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);  
    TFs_Multitaper_LOW_nK(j).values = ft_freqanalysis(cfg,datanK_bipolar);
    TFs_Multitaper_LOW_nK(j).values.freq = foi;

end



%% Do Multitaper for F

foi = 35:2.5:150; 
tw  = 0.4*ones(length(foi),1)';
fw  = 10*ones(length(foi),1)'; 
tap = floor(2.*tw.*fw);
k=2.*tw.*fw;
if floor(k)==k
    tap=k-1;
else
    tap=floor(k);
end

fRayleigh = 1./tw; 
  

for j = 1:length(subjects)
   load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.method      = 'mtmconvol';
    cfg.foi        = foi;
    cfg.t_ftimwin = tw; 
    cfg.tapsmofrq = fw;
    cfg.toi        =  -1.5:0.01:3
    cfg.pad        =    'maxperlen'; 
    cfg.taper   = 'dpss'
    
TFs_Multitaper_HIGH_eMiss(j).values = ft_freqanalysis(cfg,dataeMiss_bipolar);
TFs_Multitaper_HIGH_eMiss(j).values.freq = foi

end



for j = 1:length(subjects)
     load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
TFs_Multitaper_HIGH_nMiss(j).values = ft_freqanalysis(cfg,datanMiss_bipolar);
TFs_Multitaper_HIGH_nMiss(j).values.freq = foi


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
cfg.toi          = -1.5:0.01:3;
cfg.pad = 'maxperlen';
   
   
for j= 1:length(subjects)
     load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
    
        TFs_Multitaper_LOW_eMiss(j).values = ft_freqanalysis(cfg,dataeMiss_bipolar);
        TFs_Multitaper_LOW_eMiss(j).values.freq = foi;
      
end


for j= 1:length(subjects)
     load ([behavlist{j}, '_CleanTrialsHippo_Memory_bipolarL.mat']);
TFs_Multitaper_LOW_nMiss(j).values = ft_freqanalysis(cfg,datanMiss_bipolar);
TFs_Multitaper_LOW_nMiss(j).values.freq = foi;

end


%% concatanate high low freq remembered

 TFs_Multitaper_eR = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_eR(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_eR(j).values, TFs_Multitaper_HIGH_eR(j).values)

 end
 
 
  TFs_Multitaper_nR = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_nR(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_nR(j).values, TFs_Multitaper_HIGH_nR(j).values)
 end


%% concatanate high low freq familiar

 TFs_Multitaper_eK = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_eK(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_eK(j).values, TFs_Multitaper_HIGH_eK(j).values)

 end
 
  
  TFs_Multitaper_nK = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_nK(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_nK(j).values, TFs_Multitaper_HIGH_nK(j).values)
 end



%% concatanate high low freq forgotten

 TFs_Multitaper_eMiss = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_eMiss(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_eMiss(j).values, TFs_Multitaper_HIGH_eMiss(j).values)

 end
 
 
 
  TFs_Multitaper_nMiss = [];
 for j= 1: length(subjects)
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.appenddim = 'freq';
        TFs_Multitaper_nMiss(j).values = ft_appendfreq(cfg, TFs_Multitaper_LOW_nMiss(j).values, TFs_Multitaper_HIGH_nMiss(j).values)
 end


%% Do Baseline 

TFs_Multitaper_baseline_eR = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_eR(j).values = ft_freqbaseline(cfg,TFs_Multitaper_eR(j).values);        
end



TFs_Multitaper_baseline_nR = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_nR(j).values = ft_freqbaseline(cfg,TFs_Multitaper_nR(j).values);
end


TFs_Multitaper_baseline_eK = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_eK(j).values = ft_freqbaseline(cfg,TFs_Multitaper_eK(j).values);        
end

TFs_Multitaper_baseline_nK = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_nK(j).values = ft_freqbaseline(cfg,TFs_Multitaper_nK(j).values);        
end


TFs_Multitaper_baseline_eMiss = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_eMiss(j).values = ft_freqbaseline(cfg,TFs_Multitaper_eMiss(j).values);        
end


TFs_Multitaper_baseline_nMiss = [];
for j = 1 : length(subjects)
        cfg = [];
        cfg.baseline = [-1 -.1];
        cfg.baselinetype = 'relchange';
        TFs_Multitaper_baseline_nMiss(j).values = ft_freqbaseline(cfg,TFs_Multitaper_nMiss(j).values);
end

%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_eR(j).values.rpt= [];
    
    rj_er_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_eR(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_eR(j) = TFs_Multitaper_baseline_eR(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_eR(j).powspctrm = TFs_Multitaper_baseline_eR(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_eR(j).label = TFs_Multitaper_baseline_eR(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eR(j))
        title(['HIGH&LOW frequencies eR amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_er_hc(j).values = [rj_er_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end

%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_nR(j).values.rpt= [];
    
    rj_nr_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_nR(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_nR(j) = TFs_Multitaper_baseline_nR(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_nR(j).powspctrm = TFs_Multitaper_baseline_nR(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_nR(j).label = TFs_Multitaper_baseline_nR(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nR(j))
        title(['HIGH&LOW frequencies nR amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_nr_hc(j).values = [rj_nr_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end

%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_eK(j).values.rpt= [];
    
    rj_ek_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_eK(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_eK(j) = TFs_Multitaper_baseline_eK(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_eK(j).powspctrm = TFs_Multitaper_baseline_eK(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_eK(j).label = TFs_Multitaper_baseline_eK(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eK(j))
        title(['HIGH&LOW frequencies eK amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_ek_hc(j).values = [rj_ek_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end

%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_nK(j).values.rpt= [];
    
    rj_nk_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_nK(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_nK(j) = TFs_Multitaper_baseline_nK(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_nK(j).powspctrm = TFs_Multitaper_baseline_nK(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_nK(j).label = TFs_Multitaper_baseline_nK(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nK(j))
        title(['HIGH&LOW frequencies nK amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_nk_hc(j).values = [rj_nk_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end
%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_eK(j).values.rpt= [];
    
    rj_ek_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_eK(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_eK(j) = TFs_Multitaper_baseline_eK(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_eK(j).powspctrm = TFs_Multitaper_baseline_eK(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_eK(j).label = TFs_Multitaper_baseline_eK(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eK(j))
        title(['HIGH&LOW frequencies eK amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_ek_hc(j).values = [rj_ek_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end

%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_nK(j).values.rpt= [];
    
    rj_nk_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_nK(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_nK(j) = TFs_Multitaper_baseline_nK(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_nK(j).powspctrm = TFs_Multitaper_baseline_nK(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_nK(j).label = TFs_Multitaper_baseline_nK(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nK(j))
        title(['HIGH&LOW frequencies nK amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_nk_hc(j).values = [rj_nk_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end
%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_eMiss(j).values.rpt= [];
    
    rj_ef_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_eMiss(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j) = TFs_Multitaper_baseline_eMiss(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j).powspctrm = TFs_Multitaper_baseline_eMiss(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j).label = TFs_Multitaper_baseline_eMiss(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_eMiss(j))
        title(['HIGH&LOW frequencies eK amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_ef_hc(j).values = [rj_ef_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end

%% trial by trial visual artifact rejection in time frequency domain WITH Baseline correction 

cfg = [];
cfg.xlim =[-.5 2];
cfg.ylim = [0 150];
cfg.zlim = [-1 3];
cfg.colorbar = 'yes';

j=1;
i=1;
l=1;


for j=1:length(subjects)
    
    TFs_Multitaper_baseline_nMiss(j).values.rpt= [];
    
    rj_nk_hc(j).values = [];
    
    
    for l=1:size(TFs_Multitaper_baseline_nMiss(j).values.powspctrm,1);
        
        h = figure;
        TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j) = TFs_Multitaper_baseline_nMiss(j).values;
        TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j).powspctrm = TFs_Multitaper_baseline_nMiss(j).values.powspctrm(l,:,:,:);
        TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j).label = TFs_Multitaper_baseline_nMiss(j).values.label;
        
        
        ft_singleplotTFR(cfg,TFs_Multitaper_baseline_onechannel_onetrial_nMiss(j))
        title(['HIGH&LOW frequencies nF amy, Subject:','\_', num2str(subjects(j))])
        
        close all
            
        fprintf(strcat('******************************************trial:_',num2str(l)))
        fprintf('\n')
        artifact = 'a'; %initialization different to y or n to enter the loop
        while (artifact ~= 'y' && artifact ~= 'n')
            artifact = input('Is this trial an artifact? (y/n): ','s');
        end
        if artifact == 'y'
            
            rj_nf_hc(j).values = [rj_nf_hc(j).values l];
            
        elseif artifact == 'n'
            %do nothing
        else
            return
        end
        close all
    end
end

%% Artifact rejection in time frequency domain WITH Baseline correction 

trials2reject = [];rpt_clean=[]; 

for j=1:length(subjects)
    trials2reject.values = rj_er_hc(j).values
   
TFs_Multitaper_baseline_eR(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end  

%%
trials2reject = [];rpt_clean=[]; 

for j=1:length(subjects)
    trials2reject.values = rj_nr_hc(j).values
   
TFs_Multitaper_baseline_nR(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end  
%%

trials2reject = [];rpt_clean=[]; 

for j=1:length(subjects)
    trials2reject.values = rj_ek_hc(j).values
   
TFs_Multitaper_baseline_eK(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end  
%%  


trials2reject = [];rpt_clean=[]; 

for j=1:length(subjects)
    trials2reject.values = rj_nk_hc(j).values
   
TFs_Multitaper_baseline_nK(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end 
%% 
trials2reject = [];rpt_clean=[]; 

for j=1:length(subjects)
    trials2reject.values = rj_ef_hc(j).values
   
TFs_Multitaper_baseline_eMiss(j).values.powspctrm(trials2reject.values, :,:,:) = [];

end  
%% 
trials2reject = [];rpt_clean=[]; %re = [];rl = [];

for j=1:length(subjects)
    trials2reject.values = rj_nf_hc(j).values
   
TFs_Multitaper_baseline_nMiss(j).values.powspctrm(trials2reject.values, :,:,:) = [];

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

% save(TimeFreq_cond_9sj.mat, 'eR','eK','eF','nR','nK','nF','eKF','nKF');