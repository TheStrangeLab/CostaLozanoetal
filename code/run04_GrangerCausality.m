% This code reproduce the granger causality analysis in the direction
% amygdala to hippocampus and hippocampus to amygdala. The code reproduces
% Fig 2 a,b,c and Supplementary Fig.10

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
%%
for subj = 1:9
  cfg= [];
  dataemo{subj} = ft_appenddata(cfg, dataeR_AH{subj}, dataeK_AH{subj}, dataeF_AH{subj});
  dataneu{subj} = ft_appenddata(cfg, datanR_AH{subj}, datanK_AH{subj}, datanF_AH{subj});
  datarem{subj} = ft_appenddata(cfg, dataeR_AH{subj}, datanR_AH{subj});
  dataKforg{subj} = ft_appenddata(cfg, dataeK_AH{subj}, dataeF_AH{subj}, datanK_AH{subj}, datanF_AH{subj});
  dataeKF{subj} = ft_appenddata(cfg, dataeK_AH{subj}, dataeF_AH{subj});
  datanKF{subj} = ft_appenddata(cfg, datanK_AH{subj}, datanF_AH{subj});
end
%%

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
%%
for subj = 1:9;
    
    eR_resamp{subj}= dataeR_AH{subj};
    nR_resamp{subj}= datanR_AH{subj};
    eKF_resamp{subj}= dataeKF{subj};
    nKF_resamp{subj}= datanKF{subj}
    
    emo_resamp{subj}= dataemo{subj};
    neu_resamp{subj}= dataneu{subj};
    rem_resamp{subj}= datarem{subj};
    KF_resamp{subj}= dataKforg{subj};
    
    %% substract ensemble mean
    data_emotional_Rem{subj}=eR_resamp{subj};
    data_neutral_Rem{subj}=nR_resamp{subj};
    data_emotional_KForg{subj}=eKF_resamp{subj};
    data_neutral_KForg{subj}=nKF_resamp{subj};
    
    data_emotional{subj}=emo_resamp{subj};
    data_neutral{subj}=neu_resamp{subj};
    data_rem{subj}=rem_resamp{subj};
    data_kforg{subj}=KF_resamp{subj};
    
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
     
  
     for trial=1:size(emo_resamp{subj}.trial,2)
        dataemo_alltrials{subj}(trial,:,:)=emo_resamp{subj}.trial{trial};
    end
    
    for trial=1:size(neu_resamp{subj}.trial,2)
        dataneu_alltrials{subj}(trial,:,:)=neu_resamp{subj}.trial{trial};
    end
    
    for trial=1:size(rem_resamp{subj}.trial,2)
        datarem_alltrials{subj}(trial,:,:)=rem_resamp{subj}.trial{trial};
    end
    
    for trial=1:size(KF_resamp{subj}.trial,2)
        datakf_alltrials{subj}(trial,:,:)=KF_resamp{subj}.trial{trial};
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
    
    
  for trial=1:size(emo_resamp{subj}.trial,2)
        for chan=1:size(emo_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_emotional{subj}.trial{trial}(chan,timepoint)=emo_resamp{subj}.trial{trial}(chan,timepoint)-mean(dataemo_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    
    for trial=1:size(neu_resamp{subj}.trial,2)
        for chan=1:size(neu_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_neutral{subj}.trial{trial}(chan,timepoint)=neu_resamp{subj}.trial{trial}(chan,timepoint)-mean(dataneu_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    for trial=1:size(rem_resamp{subj}.trial,2)
        for chan=1:size(rem_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_rem{subj}.trial{trial}(chan,timepoint)=rem_resamp{subj}.trial{trial}(chan,timepoint)-mean(datarem_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
    
    for trial=1:size(KF_resamp{subj}.trial,2)
        for chan=1:size(KF_resamp{subj}.trial{trial},1)
            for timepoint=time_beginning:time_end
                data_kforg{subj}.trial{trial}(chan,timepoint)=KF_resamp{subj}.trial{trial}(chan,timepoint)-mean(datakf_alltrials{subj}(:,chan,timepoint),1);
            end
        end
    end
%%
    
    %II.III Preprocessing: 
    %-----------------------------------------------------------------------
    cfg = [];
    cfg.dftfilter ='yes';
    cfg.lpfilter = 'yes'
    cfg.lpfreq = 85
    preprocessing_emotional_Rem{subj} = ft_preprocessing(cfg,data_emotional_Rem{subj});
    preprocessing_emotional_KForg{subj}= ft_preprocessing(cfg,data_emotional_KForg{subj});
    preprocessing_neutral_Rem{subj} = ft_preprocessing(cfg,data_neutral_Rem{subj});
    preprocessing_neutral_KForg{subj}= ft_preprocessing(cfg,data_neutral_KForg{subj});
    
    preprocessing_emotional{subj} = ft_preprocessing(cfg,data_emotional{subj});
    preprocessing_neutral{subj} = ft_preprocessing(cfg,data_neutral{subj});
    preprocessing_rem{subj} = ft_preprocessing(cfg,data_rem{subj});
    preprocessing_kforg{subj} = ft_preprocessing(cfg,data_kforg{subj});
    
    % resample
    cfg            = [];
    cfg.resamplefs = 250;
    cfg.detrend    = 'no';
    preprocessing_emotional_Rem_ds{subj}  = ft_resampledata(cfg,preprocessing_emotional_Rem{subj});
    preprocessing_emotional_KForg_ds{subj}  = ft_resampledata(cfg,preprocessing_emotional_KForg{subj});
    preprocessing_neutral_Rem_ds{subj}  = ft_resampledata(cfg,preprocessing_neutral_Rem{subj});
    preprocessing_neutral_KForg_ds{subj}  = ft_resampledata(cfg,preprocessing_neutral_KForg{subj});
    
    preprocessing_emotional_ds{subj} = ft_resampledata(cfg,preprocessing_emotional{subj});
    preprocessing_neutral_ds{subj} = ft_resampledata(cfg,preprocessing_neutral{subj});
    preprocessing_rem_ds{subj} = ft_resampledata(cfg,preprocessing_rem{subj});
    preprocessing_kforg_ds{subj} = ft_resampledata(cfg,preprocessing_kforg{subj});
    
    
    cfg         = [];
    cfg.order   = 9;
    cfg.toolbox = 'bsmart';
    cfg.toi       = -0.5:0.01:1.5;
    cfg.t_ftimwin = 0.4;
    %interaction
    mdata_emotional_Rem{subj} = ft_mvaranalysis(cfg, preprocessing_emotional_Rem_ds{subj});
    mdata_emotional_KForg{subj} = ft_mvaranalysis(cfg, preprocessing_emotional_KForg_ds{subj});
    mdata_neutral_Rem{subj} = ft_mvaranalysis(cfg, preprocessing_neutral_Rem_ds{subj});
    mdata_neutral_KForg{subj} = ft_mvaranalysis(cfg, preprocessing_neutral_KForg_ds{subj});
    
    %main effects
    mdata_emotional{subj} = ft_mvaranalysis(cfg, preprocessing_emotional_ds{subj});
    mdata_neutral{subj} = ft_mvaranalysis(cfg, preprocessing_neutral_ds{subj});
    mdata_rem{subj} = ft_mvaranalysis(cfg, preprocessing_rem_ds{subj});
    mdata_kforg{subj} = ft_mvaranalysis(cfg, preprocessing_kforg_ds{subj});
    
        
    %     I.III Analysis with convolution method
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    freq_choice=2:80;
    
    
    cfg = [];
    cfg.method = 'mvar';
    cfg.foi        = freq_choice;
    %interaction
    Freq_emotional_Rem{subj} =  ft_freqanalysis_mvar(cfg,mdata_emotional_Rem{subj});
    Freq_emotional_KForg{subj} =  ft_freqanalysis_mvar(cfg,mdata_emotional_KForg{subj});
    Freq_neutral_Rem{subj} =  ft_freqanalysis_mvar(cfg,mdata_neutral_Rem{subj});
    Freq_neutral_KForg{subj} =  ft_freqanalysis_mvar(cfg,mdata_neutral_KForg{subj});
    % main effects
    Freq_emotional{subj} =  ft_freqanalysis_mvar(cfg,mdata_emotional{subj} );
    Freq_neutral{subj} =  ft_freqanalysis_mvar(cfg,mdata_neutral{subj});
    Freq_rem{subj} =  ft_freqanalysis_mvar(cfg,mdata_rem{subj});
    Freq_kforg{subj} =  ft_freqanalysis_mvar(cfg,mdata_kforg{subj});
    
    
    cfg           = [];
    cfg.method    = 'granger';
    % interaction
    emotional_Rem_granger{subj}  = ft_connectivityanalysis(cfg, Freq_emotional_Rem{subj});
    emotional_KForg_granger{subj}  = ft_connectivityanalysis(cfg,Freq_emotional_KForg{subj});
    neutral_Rem_granger{subj}  = ft_connectivityanalysis(cfg, Freq_neutral_Rem{subj});
    neutral_KForg_granger{subj}  = ft_connectivityanalysis(cfg,Freq_neutral_KForg{subj});
    
    %main effects
    emotional_granger{subj}  = ft_connectivityanalysis(cfg, Freq_emotional{subj});
    neutral_granger{subj}  = ft_connectivityanalysis(cfg,Freq_neutral{subj});
    rem_granger{subj}  = ft_connectivityanalysis(cfg, Freq_rem{subj});
    kforg_granger{subj}  = ft_connectivityanalysis(cfg,Freq_kforg{subj});
    
    % interaction
    GA_emotional_Rem_HCtoAmyg(subj,:,:)= emotional_Rem_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_emotional_KForg_HCtoAmyg(subj,:,:)=emotional_KForg_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_emotional_Rem_AmygtoHC(subj,:,:)=emotional_Rem_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_emotional_KForg_AmygtoHC(subj,:,:)=emotional_KForg_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
    GA_neutral_Rem_HCtoAmyg(subj,:,:)=neutral_Rem_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_neutral_KForg_HCtoAmyg(subj,:,:)=neutral_KForg_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_neutral_Rem_AmygtoHC(subj,:,:)=neutral_Rem_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_neutral_KForg_AmygtoHC(subj,:,:)=neutral_KForg_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
    % main effects
    GA_emotional_HCtoAmyg(subj,:,:)= emotional_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_neutral_HCtoAmyg(subj,:,:)=neutral_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_emotional_AmygtoHC(subj,:,:)=emotional_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_neutral_AmygtoHC(subj,:,:)=neutral_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    
    GA_rem_HCtoAmyg(subj,:,:)=rem_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_kforg_HCtoAmyg(subj,:,:)=kforg_granger{subj}.grangerspctrm(HC_list{subj},Amy_list{subj},:,:);
    GA_rem_AmygtoHC(subj,:,:)=rem_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);
    GA_kforg_AmygtoHC(subj,:,:)=kforg_granger{subj}.grangerspctrm(Amy_list{subj},HC_list{subj},:,:);        

end

%% figure
clear all
close all
clc

subjects=[60 13 15 16 160 25 32 10 33];
mlist={'s60', 's13', 's15', 's16', 's160', 's25', 's32', 's10', 's33'};
behavlist = {'Patient60','Patient13','Patient15','Patient16','Patient160','Patient25','Patient32','PatientZ_10','Patient33'};

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','Coherence','GA_LFPs.mat'))

% cd C:\Users\manuela\Desktop\data2github\Coherence
% load GA_LFPs.mat
% 
% cd cd C:\Users\manuela\Desktop\data2github\Granger
% load Granger_Var_timefreq_Int.mat
% load Freq_mvar_Int.mat
% load Granger_Var_timefreq_meffects.mat
% load Freq_mvar_meffects.mat
%%
GTFs_Multitaper_baseline_Unpl.time = []
GTFs_Multitaper_baseline_Unpl.time = Freq_rem{1}.time
GTFs_Multitaper_baseline_Unpl.freq = []
GTFs_Multitaper_baseline_Unpl.freq = Freq_rem{1}.freq
GTFs_Multitaper_baseline_Unpl.label = []
GTFs_Multitaper_baseline_Unpl.label = {'AH'}


GA_stat_emotional_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_emotional_HCtoAmyg.powspctrm=[];
GA_stat_emotional_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_emotional_HCtoAmyg(:,:,:,:),4));

GA_stat_neutral_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_neutral_HCtoAmyg.powspctrm=[];
GA_stat_neutral_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_neutral_HCtoAmyg(:,:,:,:),4));

GA_stat_rem_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_rem_HCtoAmyg.powspctrm=[];
GA_stat_rem_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_rem_HCtoAmyg(:,:,:,:),4));

GA_stat_kforg_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_kforg_HCtoAmyg.powspctrm=[];
GA_stat_kforg_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_kforg_HCtoAmyg(:,:,:,:),4));
%%
GA_stat_emotional_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_emotional_AmygtoHC.powspctrm=[];
GA_stat_emotional_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_emotional_AmygtoHC(:,:,:,:),4));

GA_stat_neutral_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_neutral_AmygtoHC.powspctrm=[];
GA_stat_neutral_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_neutral_AmygtoHC(:,:,:,:),4));

GA_stat_rem_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_rem_AmygtoHC.powspctrm=[];
GA_stat_rem_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_rem_AmygtoHC(:,:,:,:),4));

GA_stat_kforg_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_kforg_AmygtoHC.powspctrm=[];
GA_stat_kforg_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_kforg_AmygtoHC(:,:,:,:),4));

%%
GA_stat_direction_Emo= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_Emo.powspctrm =[];
GA_stat_direction_Emo.powspctrm(:,1,:,:)=GA_stat_emotional_AmygtoHC.powspctrm - GA_stat_emotional_HCtoAmyg.powspctrm;

GA_stat_direction_Neu= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_Neu.powspctrm =[];
GA_stat_direction_Neu.powspctrm(:,1,:,:)=GA_stat_neutral_AmygtoHC.powspctrm - GA_stat_neutral_HCtoAmyg.powspctrm;

GA_stat_direction_Rem= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_Rem.powspctrm =[];
GA_stat_direction_Rem.powspctrm(:,1,:,:)=GA_stat_rem_AmygtoHC.powspctrm - GA_stat_rem_HCtoAmyg.powspctrm;

GA_stat_direction_KForg= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_KForg.powspctrm =[];
GA_stat_direction_KForg.powspctrm(:,1,:,:)=GA_stat_kforg_AmygtoHC.powspctrm - GA_stat_kforg_HCtoAmyg.powspctrm;

%%
cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];%
cfg.ylim = [2 35];%'maxmin';
cfg.zlim = [-0.015 0.015];%[0 0.2];%'maxmin';
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,GA_stat_direction_Emo);
tit=title('Unpleasant (A2H-H2A)');set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold');
hold on
subplot(2,2,2) ; ft_singleplotTFR(cfg,GA_stat_direction_Neu);
tit=title('Neutral (A2H-H2A)');set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold');

%%
cfg = [];
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];%
cfg.ylim = [2 35];%'maxmin';
cfg.zlim = [-0.015 0.015];%[0 0.2];%'maxmin';
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,GA_stat_direction_Rem);
tit=title('Rem (A2H-H2A)');set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold');
hold on
subplot(2,2,2) ; ft_singleplotTFR(cfg,GA_stat_direction_KForg);
tit=title('Kforg (A2H-H2A)');set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold');

%%
GA_stat_emoa2h = GA_stat_emotional_AmygtoHC
GA_stat_emoa2h.powspctrm = []
GA_stat_emoa2h.powspctrm = GA_stat_emotional_AmygtoHC.powspctrm - GA_stat_neutral_AmygtoHC.powspctrm

GA_stat_mema2h = GA_stat_rem_AmygtoHC
GA_stat_mema2h.powspctrm = []
GA_stat_mema2h.powspctrm = GA_stat_rem_AmygtoHC.powspctrm - GA_stat_kforg_AmygtoHC.powspctrm
%%
cfg = []; 
cfg.figure = 'gcf'
cfg.xlim = [-.5 1.5];%
cfg.ylim = [2 35];%'maxmin';
cfg.zlim = [-0.015 0.015];%[0 0.2];%'maxmin';
h = figure;set(h, 'Position', get(0, 'Screensize'))
subplot(2,2,1) ; ft_singleplotTFR(cfg,GA_stat_emoa2h);
tit=title('Unpleasant - Neutral A2H');set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold');
hold on
subplot(2,2,2) ; ft_singleplotTFR(cfg,GA_stat_mema2h);
tit=title('Remember - Kforgotten A2H');set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold');


%% cluster stat

cfg=[];
cfg.method = 'montecarlo'%
cfg.statistic = 'depsamplesT';

cfg.correctm = 'cluster';

t1=0;
t2=1;
f1=2;
f2=34;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];

cfg.tail             = 0; 
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;
cfg.numrandomization = 'all';%10000;

cfg.neighbours = [];
cfg.ivar = 1;
cfg.uvar = 2;

cfg.design = [ones(1,9) ones(1,9).*2;[1:9] [1:9]];

emo_A2H = ft_freqstatistics(cfg, GA_stat_emotional_AmygtoHC, GA_stat_neutral_AmygtoHC);
emo_H2A = ft_freqstatistics(cfg, GA_stat_emotional_HCtoAmyg, GA_stat_neutral_HCtoAmyg);

mem_A2H = ft_freqstatistics(cfg, GA_stat_rem_AmygtoHC, GA_stat_kforg_AmygtoHC);
mem_H2A = ft_freqstatistics(cfg, GA_stat_rem_HCtoAmyg, GA_stat_kforg_HCtoAmyg);

emodir = ft_freqstatistics(cfg, GA_stat_direction_Emo, GA_stat_direction_Neu);
memdir = ft_freqstatistics(cfg, GA_stat_direction_Rem, GA_stat_direction_KForg);

%% 
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1]
cfg.ylim = [2 35]
cfg.zlim = [-4 4];

E = emo_A2H;
E.powspctrm = emo_A2H.stat.*emo_A2H.mask;
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('emo A2H'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title(['a2h  p=' num2str([emo_A2H.posclusters(1).prob])])

%% This code reproduce Fig.2a
h=figure;set(h, 'Position', get(0, 'Screensize'))

logRelative1 = emo_A2H
logRelative1.mask = emo_A2H.mask;
logRelative1.powspctrm = emo_A2H.stat
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
ft_singleplotTFR(cfg,logRelative1);

%% Figure 2c
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1]
cfg.ylim = [2 35]
cfg.zlim = [-4 4];

E = emo_H2A;
E.powspctrm = emo_H2A.stat
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('emo H2A'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title(['h2a  p=' num2str([emo_H2A.posclusters(1).prob])])

%% supplementary Fig 10c
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1]
cfg.ylim = [2 35]
cfg.zlim = [-4 4];

E = mem_H2A;
E.powspctrm = mem_H2A.stat
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('mem h2a'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title(['h2a  p=' num2str([mem_H2A.posclusters(1).prob])])
%% supplementary Fig 10a
h=figure;set(h, 'Position', get(0, 'Screensize'))
cfg = [];
cfg.xlim =[0 1]
cfg.ylim = [2 35]
cfg.zlim = [-4 4];

E = mem_A2H;
E.powspctrm = mem_A2H.stat
E.dimord = 'chan_freq_time';
ft_singleplotTFR(cfg,E); colorbar;  tit=title('mem ah2'); xlabel('time in s'); ylabel('frequency in Hz');
set(findobj(tit,'type','text'),'fontsize',8,'FontWeight','bold') ; 
title(['a2h  p=' num2str(min([mem_A2H.posclusters(1).prob,mem_A2H.negclusters(1).prob]))])

%%
maskt= sum(squeeze(emo_A2H.mask),1);
tvec = emo_A2H.time(maskt>0);
toi=[min(tvec) max(tvec)]

maskf= sum(squeeze(emo_A2H.mask),2);
fvec = emo_A2H.freq(maskf>0);
foi=[min(fvec) max(fvec)]

t=toi
f=foi

pt1 = nearest(GA_stat_emotional_AmygtoHC.time,t(1));
pt2 = nearest(GA_stat_emotional_AmygtoHC.time,t(2));
pf1 = nearest(GA_stat_emotional_AmygtoHC.freq,f(1));
pf2 = nearest(GA_stat_emotional_AmygtoHC.freq,f(2));

% create matrix of conditions 
mat(:,1) = squeeze(mean(mean(GA_stat_emotional_AmygtoHC.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));
mat(:,2) = squeeze(mean(mean(GA_stat_neutral_AmygtoHC.powspctrm(:,:,pf1:pf2,pt1:pt2),4),3));

%%
oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'))

x = [ones(9,1) ones(9,1)*2];
y = [mat(:,1) mat(:,2)];
figure;beeswarm(x(:),y(:),'sort_style','up','dot_size',4,'overlay_style','sd','colormap',[1 0 0; 0 0 0])
ylim([-0.004 0.038]);
xticklabels({[],'emotional a2h',[],'neutral a2h'})
grid off
hold on

%% calculate the interaction
GTFs_Multitaper_baseline_Unpl.time = []
GTFs_Multitaper_baseline_Unpl.time = Freq_neutral_Rem{1}.time
GTFs_Multitaper_baseline_Unpl.freq = []
GTFs_Multitaper_baseline_Unpl.freq = Freq_neutral_Rem{1}.freq
GTFs_Multitaper_baseline_Unpl.label = []
GTFs_Multitaper_baseline_Unpl.label = {'AH'}


GA_stat_emotional_Rem_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_emotional_Rem_HCtoAmyg.powspctrm=[];
GA_stat_emotional_Rem_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_emotional_Rem_HCtoAmyg([1:7,9],:,:,:),4));

GA_stat_emotional_KForg_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_emotional_KForg_HCtoAmyg.powspctrm=[];
GA_stat_emotional_KForg_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_emotional_KForg_HCtoAmyg([1:7,9],:,:,:),4));

GA_stat_neutral_Rem_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_neutral_Rem_HCtoAmyg.powspctrm=[];
GA_stat_neutral_Rem_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_neutral_Rem_HCtoAmyg([1:7,9],:,:,:),4));

GA_stat_neutral_KForg_HCtoAmyg=GTFs_Multitaper_baseline_Unpl;
GA_stat_neutral_KForg_HCtoAmyg.powspctrm=[];
GA_stat_neutral_KForg_HCtoAmyg.powspctrm(:,1,:,:)=squeeze(mean(GA_neutral_KForg_HCtoAmyg([1:7,9],:,:,:),4));
%%
GA_stat_emotional_Rem_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_emotional_Rem_AmygtoHC.powspctrm=[];
GA_stat_emotional_Rem_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_emotional_Rem_AmygtoHC([1:7,9],:,:,:),4));

GA_stat_emotional_KForg_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_emotional_KForg_AmygtoHC.powspctrm=[];
GA_stat_emotional_KForg_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_emotional_KForg_AmygtoHC([1:7,9],:,:,:),4));

GA_stat_neutral_Rem_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_neutral_Rem_AmygtoHC.powspctrm=[];
GA_stat_neutral_Rem_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_neutral_Rem_AmygtoHC([1:7,9],:,:,:),4));

GA_stat_neutral_KForg_AmygtoHC=GTFs_Multitaper_baseline_Unpl;
GA_stat_neutral_KForg_AmygtoHC.powspctrm=[];
GA_stat_neutral_KForg_AmygtoHC.powspctrm(:,1,:,:)=squeeze(mean(GA_neutral_KForg_AmygtoHC([1:7,9],:,:,:),4));


%% for interaction effects

GA_stat_direction_emotional_KForg= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_emotional_KForg.powspctrm =[];
GA_stat_direction_emotional_KForg.powspctrm(:,1,:,:)=GA_stat_emotional_KForg_AmygtoHC.powspctrm - GA_stat_emotional_KForg_HCtoAmyg.powspctrm;

GA_stat_direction_emotional_Rem= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_emotional_Rem.powspctrm =[];
GA_stat_direction_emotional_Rem.powspctrm(:,1,:,:)=GA_stat_emotional_Rem_AmygtoHC.powspctrm - GA_stat_emotional_Rem_HCtoAmyg.powspctrm;

GA_stat_direction_neutral_KForg= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_neutral_KForg.powspctrm =[];
GA_stat_direction_neutral_KForg.powspctrm(:,1,:,:)=GA_stat_neutral_KForg_AmygtoHC.powspctrm - GA_stat_neutral_KForg_HCtoAmyg.powspctrm;

GA_stat_direction_neutral_Rem= GTFs_Multitaper_baseline_Unpl;
GA_stat_direction_neutral_Rem.powspctrm =[];
GA_stat_direction_neutral_Rem.powspctrm(:,1,:,:)=GA_stat_neutral_Rem_AmygtoHC.powspctrm - GA_stat_neutral_Rem_HCtoAmyg.powspctrm;
%%
URUKF = GTFs_Multitaper_baseline_Unpl;
URUKF.powspctrm = []
URUKF.powspctrm =GA_stat_direction_emotional_Rem.powspctrm - GA_stat_direction_emotional_KForg.powspctrm

NRNKF = GTFs_Multitaper_baseline_Unpl;
NRNKF.powspctrm = []
NRNKF.powspctrm =GA_stat_direction_neutral_Rem.powspctrm - GA_stat_direction_neutral_KForg.powspctrm
%%
URUKF_a2h = GTFs_Multitaper_baseline_Unpl;
URUKF_a2h.powspctrm = []
URUKF_a2h.powspctrm =GA_stat_emotional_Rem_AmygtoHC.powspctrm - GA_stat_emotional_KForg_AmygtoHC.powspctrm;

NRNKF_a2h = GTFs_Multitaper_baseline_Unpl;
NRNKF_a2h.powspctrm = []
NRNKF_a2h.powspctrm =GA_stat_neutral_Rem_AmygtoHC.powspctrm - GA_stat_neutral_KForg_AmygtoHC.powspctrm;
%%
URUKF_h2a = GTFs_Multitaper_baseline_Unpl;
URUKF_h2a.powspctrm = []
URUKF_h2a.powspctrm =GA_stat_emotional_Rem_HCtoAmyg.powspctrm - GA_stat_emotional_KForg_HCtoAmyg.powspctrm;

NRNKF_h2a = GTFs_Multitaper_baseline_Unpl;
NRNKF_h2a.powspctrm = []
NRNKF_h2a.powspctrm =GA_stat_neutral_Rem_HCtoAmyg.powspctrm - GA_stat_neutral_KForg_HCtoAmyg.powspctrm;
%% cluter stat

cfg=[];
cfg.method = 'montecarlo'
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';

t1=0;
t2=1;
f1=2;
f2=34;
cfg.latency = [t1 t2];
cfg.frequency = [f1 f2];


cfg.tail             = 0; 
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusteralpha     = 0.05;%
cfg.numrandomization = 'all';%10000;

cfg.neighbours = []; 
cfg.ivar = 1;
cfg.uvar = 2;

cfg.design = [ones(1,8) ones(1,8).*2;[1:8] [1:8]];
int_dir = ft_freqstatistics(cfg, URUKF, NRNKF);
int_a2h = ft_freqstatistics(cfg, URUKF_a2h, NRNKF_a2h);
int_h2a = ft_freqstatistics(cfg, URUKF_h2a, NRNKF_h2a);
