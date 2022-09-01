clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

subjects=[60 13 15 16 160 25 32 10 33];
mlist={'s60', 's13', 's15', 's16', 's160', 's25', 's32', 's10', 's33'};
behavlist = {'Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient32','PatientZ_10','Patient33'};

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','DataCleanAH_RKF.mat'))

%s60
Amy_list{1}=[2];
HC_list{1}=[3];
%s13
Amy_list{2}=[2];
HC_list{2}=[4];
%s15
Amy_list{3}=[2];
HC_list{3}=[3];
%s16
Amy_list{4}=[1];
HC_list{4}=[3];
%s160
Amy_list{5}=[1];
HC_list{5}=[3];
%s25
Amy_list{6}=[1];
HC_list{6}=[2];
%s32
Amy_list{7}=[1];
HC_list{7}=[2];
%s10
Amy_list{8}=[1];
HC_list{8}=[2];
%s33
Amy_list{9}=[2];
HC_list{9}=[3];

for subj = 1:9
  cfg= [];
  dataeKF{subj} = ft_appenddata(cfg, dataeK_AH{subj}, dataeF_AH{subj});
  datanKF{subj} = ft_appenddata(cfg, datanK_AH{subj}, datanF_AH{subj});
end

for subj = 1:9;
    
    eR_resamp{subj}= dataeR_AH{subj};
    nR_resamp{subj}= datanR_AH{subj};
    eKF_resamp{subj}= dataeKF{subj};
    nKF_resamp{subj}= datanKF{subj}
    %% substract ensemble mean
    data_emotional_Rem{subj}=eR_resamp{subj};
    data_neutral_Rem{subj}=nR_resamp{subj};
    data_emotional_KForg{subj}=eKF_resamp{subj};
    data_neutral_KForg{subj}=nKF_resamp{subj};

    time_beginning=1;
    time_end=size(eR_resamp{subj}.time{1},2);
    
    for trial=1:size(eR_resamp{subj}.trial,2)
        dataeR_alltrials{subj}(trial,:,:)=eR_resamp{subj}.trial{trial};
    end
    
    for trial=1:size(nR_resamp{subj}.trial,2)
        datanR_alltrials{subj}(trial,:,:)=nR_resamp{subj}.trial{trial};
    end
    
    for trial=1:size(eKF_resamp{subj}.trial,2)
        dataeKF_alltrials{subj}(trial,:,:)=eKF_resamp{subj}.trial{trial};
    end
    
     for trial=1:size(nKF_resamp{subj}.trial,2)
        datanKF_alltrials{subj}(trial,:,:)=nKF_resamp{subj}.trial{trial};
     end
    
    for trial=1:size(eR_resamp{subj}.trial,2)
        for chan=1:size(eR_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_emotional_Rem{subj}.trial{trial}(chan,timepoint)=eR_resamp{subj}.trial{trial}(chan,timepoint)-mean(dataeR_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    
    for trial=1:size(nR_resamp{subj}.trial,2)
        for chan=1:size(nR_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_neutral_Rem{subj}.trial{trial}(chan,timepoint)=nR_resamp{subj}.trial{trial}(chan,timepoint)-mean(datanR_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    for trial=1:size(eKF_resamp{subj}.trial,2)
        for chan=1:size(eKF_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_emotional_KForg{subj}.trial{trial}(chan,timepoint)=eKF_resamp{subj}.trial{trial}(chan,timepoint)-mean(dataeKF_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    for trial=1:size(nKF_resamp{subj}.trial,2)
        for chan=1:size(nKF_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_neutral_KForg{subj}.trial{trial}(chan,timepoint)=nKF_resamp{subj}.trial{trial}(chan,timepoint)-mean(datanKF_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    cfg= []
    data_emotional{subj} =  ft_appenddata(cfg, data_emotional_Rem{subj},  data_emotional_KForg{subj});
    data_neutral{subj} =  ft_appenddata(cfg, data_neutral_Rem{subj},  data_neutral_KForg{subj});
    data_remember{subj} =  ft_appenddata(cfg, data_emotional_Rem{subj},  data_neutral_Rem{subj});
    data_kforg{subj} =  ft_appenddata(cfg, data_emotional_KForg{subj},  data_neutral_KForg{subj});
    
    %II.III Preprocessing: Filtern
    %-----------------------------------------------------------------------
    cfg = [];
    cfg.dftfilter ='yes';
    preprocessing_emotional_Rem{subj} = ft_preprocessing(cfg,data_emotional_Rem{subj});
    preprocessing_emotional_KForg{subj}= ft_preprocessing(cfg,data_emotional_KForg{subj});

    preprocessing_neutral_Rem{subj} = ft_preprocessing(cfg,data_neutral_Rem{subj});
    preprocessing_neutral_KForg{subj}= ft_preprocessing(cfg,data_neutral_KForg{subj});

  
    foi = 2:140; % Frequency of interest
    cfg              = [];
    cfg.output       = 'fourier';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';

    cfg.foi          = foi;
    cfg.t_ftimwin    = 4./cfg.foi; 
    cfg.toi          = -0.5:0.01:1.5;% time window "slides" in steps of 0.01 sec (10 ms)
    cfg.pad = 'maxperlen';
    
    Freq_emotional_Rem{subj} =  ft_freqanalysis(cfg,preprocessing_emotional_Rem{subj});
    Freq_emotional_KForg{subj} =  ft_freqanalysis(cfg,preprocessing_emotional_KForg{subj});

    Freq_neutral_Rem{subj} =  ft_freqanalysis(cfg,preprocessing_neutral_Rem{subj});
    Freq_neutral_KForg{subj} =  ft_freqanalysis(cfg,preprocessing_neutral_KForg{subj});

    cfg           = [];
    cfg.method    = 'coh'
    emotional_Rem_coh{subj}  = ft_connectivityanalysis(cfg, Freq_emotional_Rem{subj});
    emotional_KForg_coh{subj}  = ft_connectivityanalysis(cfg,Freq_emotional_KForg{subj});

    neutral_Rem_coh{subj}  = ft_connectivityanalysis(cfg, Freq_neutral_Rem{subj});
    neutral_KForg_coh{subj}  = ft_connectivityanalysis(cfg,Freq_neutral_KForg{subj});


   GA_emotional_Rem_AmyHC(subj,:,:)=emotional_Rem_coh{subj}.cohspctrm(Amy_list{subj},HC_list{subj},:,:);
   GA_emotional_KForg_AmyHC(subj,:,:)=emotional_KForg_coh{subj}.cohspctrm(Amy_list{subj},HC_list{subj},:,:);

   GA_neutral_Rem_AmyHC(subj,:,:)=neutral_Rem_coh{subj}.cohspctrm(Amy_list{subj},HC_list{subj},:,:);
   GA_neutral_KForg_AmyHC(subj,:,:)=neutral_KForg_coh{subj}.cohspctrm(Amy_list{subj},HC_list{subj},:,:);

     
% 
%   save('Coherence_Var_timefreq.mat', 'GA_emotional_Rem_AmyHC', 'GA_neutral_Rem_AmyHC', 'GA_emotional_KForg_AmyHC','GA_neutral_KForg_AmyHC');
%   save ('Freq_mvar.mat','Freq_emotional_Rem', 'Freq_emotional_KForg','Freq_neutral_Rem', 'Freq_neutral_KForg');
end
%% This reproduce Supplementary Fig. 9 a and b
clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20190203\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','Coherence','Coherence_Var_timefreq.mat'))
load(fullfile(oripath,'data2github','Coherence','Freq_mvar'))
load(fullfile(oripath,'data2github','Coherence','GA_LFPs.mat'))


Stat = GaLFPs
Stat.individual=[];
Stat.time=[];


Stat.time= 2:140;
%%
Nsubjects=size(GA_neutral_Rem_AmyHC(:,:,:),1);

%% main effect emo


GA_neutral=Stat;

GA_neutral.individual=[];

mean_neutral(1,[1:7,9],1,:)= mean(GA_neutral_Rem_AmyHC([1:7,9],:,:),3);

mean_neutral(2,:,1,:)= mean(GA_neutral_KForg_AmyHC(:,:,:),3);

GA_neutral.individual(:,1,:)=squeeze(mean(mean_neutral,1));



GA_emotional=Stat;

GA_emotional.individual=[];

mean_emotional(1,:,1,:)=mean(GA_emotional_Rem_AmyHC(:,:,:),3);

mean_emotional(2,:,1,:)=mean(GA_emotional_KForg_AmyHC(:,:,:),3);

GA_emotional.individual(:,1,:)=squeeze(mean(mean_emotional,1));


%% main effect memory success (outcome)


GA_Rem=Stat;

GA_Rem.individual=[];

mean_Rem(1,:,1,:)=mean(GA_emotional_Rem_AmyHC(:,:,:),3);

mean_Rem(2,[1:7,9],1,:)=mean(GA_neutral_Rem_AmyHC([1:7,9],:,:),3);

GA_Rem.individual(:,1,:)=squeeze(mean(mean_Rem,1));



GA_KForg=Stat;

GA_KForg.individual=[];

mean_KForg(1,:,1,:)=mean(GA_emotional_KForg_AmyHC(:,:,:),3);

mean_KForg(2,:,1,:)=mean(GA_neutral_KForg_AmyHC(:,:,:),3);

GA_KForg.individual(:,1,:)=squeeze(mean(mean_KForg,1));

%% Plot results 

emo = squeeze(GA_emotional.individual)
neu = squeeze(GA_neutral.individual)
rem = squeeze(GA_Rem.individual)
kforg = squeeze(GA_KForg.individual)

norm_emo = normalize(emo)
norm_neu = normalize(neu)

norm_rem = normalize(rem)
norm_kforg = normalize(kforg)
%%
for i=1:139
avgA = mean(norm_emo);
stdA(:,i) = std(avgA)/sqrt(9)
end

for i=1:139
avgB = mean(norm_neu);
stdB(:,i) = std(avgB)/sqrt(9)
end

figure
x= 1:1:139
shadedErrorBar(x,avgA,stdA,'r')
hold on
shadedErrorBar(x,avgB,stdB,'b')
hold on
legend('emo', 'neu')
xlabel('frequency hz')
ylabel ('coherence')
box off

clear avgA, clear stdA
for i=1:139
avgA = mean(norm_rem);
stdA(:,i) = std(avgA)/sqrt(9)
end

clear avgB, clear stdB

for i=1:139
avgB = mean(norm_kforg);
stdB(:,i) = std(avgB)/sqrt(9)
end

figure
x= 1:1:139
shadedErrorBar(x,avgA,stdA,'m')
hold on
shadedErrorBar(x,avgB,stdB,'k')
hold on
legend('rem', 'kforg')
xlabel('frequency hz')
ylabel ('coherence')
box off

%%
Stat.individual=[];

GA_neutral=Stat;
GA_neutral.individual=[];
GA_neutral.individual(:,1,:)= norm_neu;

GA_emotional=Stat;
GA_emotional.individual=[];
GA_emotional.individual(:,1,:)= norm_emo;
%%
GA_Rem=Stat;
GA_Rem.individual=[];
GA_Rem.individual(:,1,:)= norm_rem;

GA_KForg=Stat;
GA_KForg.individual=[];
GA_KForg.individual(:,1,:)= norm_kforg;
%%
cfg = [];

cfg.method = 'montecarlo';

cfg.statistic = 'depsamplesT';

cfg.correctm = 'cluster';

cfg.clusteralpha = 0.05; 

cfg.tail = 0; 

cfg.clustertail = 0;

cfg.alpha = 0.025;

cfg.numrandomization =10000; 


Nsubjects = 9
%Design matrix

%-------------------

design = zeros(2,2*Nsubjects);

for i = 1:Nsubjects

    design(1,i) = i;

end

for i = 1:Nsubjects

    design(1,Nsubjects+i) = i;

end

design(2,1:Nsubjects)= 1;

design(2,Nsubjects+1:2*Nsubjects) = 2;



cfg.design = design;

cfg.uvar  = 1;  

cfg.ivar  = 2;

main_effect_emotion=ft_timelockstatistics(cfg,GA_emotional,GA_neutral);
main_effect_memory=ft_timelockstatistics(cfg,GA_Rem,GA_KForg);
