clear all; close all; clc;

% add fieldtrip
restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','eeg_toolbox_v1_3_2'))% kahana eeg toolbox
addpath(fullfile(oripath,'Costalozanoetal','code','utils','MatlabImportExport_v6.0.0'))% 
warning('off', 'all'); % suppress all warnings

% select manually depending on wether you want to look at the data 
% in the positive or negative spikes folder , either cluster 1 or cluster 2

%direction='posSpikes';
direction='negSpikes';

cluster= 'Cluster1' % when you change this check iclus value either 1 or 2
%cluster= 'Cluster2' 

bins = '20_bins' 


% direct to the folder with wave clus results
spike_path      = ['C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\SpikeExtraction_20210210\' direction];


% settings
angle_bins      = deg2rad(-180:18:180);
angle_centers  	= transpose(movmean(angle_bins, 2, 'endpoints', 'discard'));
t_pre_stim = 1% How many seconds before the trigger to take - scenes were presented for 500 ms + 3500 of iti
t_post_stim = 3;% How many seconds after the trigger to take


trans_width     = 0.15;
idealresponse   = [0, 0, 1, 1, 0, 0];
numSurrogates   = 1000;
cutoffHighPower = 0;


% analysis options
qRemoveSpikes           = 'NoremSpikes'; % question whether to remove spikes from the LFP before processing it: 'remSpikes' or 'notRemSpikes'
data4Rayleigh           = 'raw4Rayleigh'; % whether to normalize data before Rayleigh test: 'raw4Rayleigh' or 'norm4Rayleigh'
surroCreationType       = 'fullyRandom'; % how to create the surrogate data: 'fullyRandom' or 'closeAround'

% subjects list
subjects    = {...
    'Patient1'; ...
    'Patient2'; ...
    'Patient200'; ...
    };

pz = [1 2 200];
%% get LFP from the amygdala
chanofinterest{1} = {'uAL1','uAL2'};
labelorg{1} = {'uAL1','uAL2'};
bipolabel{1} = {'uAL1 - uAL2'};
bipolmat{1} = [+1 -1];

chanofinterest{2} = {'uAL1','uAL2'};
labelorg{2} = {'uAL1','uAL2'};
bipolabel{2} = {'uAL1 - uAL2'};
bipolmat{2} = [+1 -1];

chanofinterest{3} = {'uAR1','uAR2'};
labelorg{3} = {'uAR1','uAR2'};
bipolabel{3} = {'uAR1 - uAR2'};
bipolmat{3} = [+1 -1];

v=0
%% loop through subjects
for isub = 1:length(subjects)
    v=v+1;
    
    fprintf('\n\nSUBJECT: %s.\n', subjects{isub});
    rng(1); % for reproducibility
    
    
    
    %% get onsets
    
    load(fullfile(oripath,'data2github','SingleUnit',subjects{isub},'Behaviour','enc_onsets.mat'))
   
    clear hdr
    clear evR
    clear evCorrected
    clear evClean
    clear trl
    clear data
    clear data_eR
    clear data_eKF
    clear ekf
    clear onsets_for_spikes
    clear onsets_for_eRspikes
    clear onsets_for_eKFspikes

    %  This part of the script process the raw data that are provided in
    %  data2github folder and can be load from line 164
    
    %   set path for the macrowires and read the data
    
%         strNLXFolderPath = ['~\0000_IAPS_Spike\',subjects{isub},'\EEGData\u\']    
%         hdr = ft_read_header(strNLXFolderPath);
%         evR = ft_read_event(strNLXFolderPath);
%         indEventRangeStart = 1; % Beginning of event structure % CHANGE THIS
%         indEventRangeStop = length(evR)-1; % End of event structure % CHANGE THIS
%         evSelected = evR(indEventRangeStart:indEventRangeStop);
%         evClean = Get_Clean_Event_Structure(evSelected);
%         strFirstFileName = dir(fullfile(strNLXFolderPath,'*.ncs'));
%         strFirstFileName = {strFirstFileName.name}';
%         strFirstFileName = strFirstFileName(1);
%         strFirstFilePath = strcat(strNLXFolderPath,strFirstFileName);
%         ts = ft_read_data(strFirstFilePath,'timestamp', true);
%         evCorrected = evClean;
%         for ii = 1:length(evCorrected)
%             % The sample corresponding to this timestamp
%             [~,indSample] = min(abs(double(ts)-double(evCorrected(ii).timestamp)));
%             % Update field 'sample'
%             evCorrected(ii).sample = indSample;
%         end
%     
%         cfg = []; % initialize
%         cfg.dataset = strNLXFolderPath;
%         cfg.trialdef.eventtype  = 'trigger'; % Define trials based on triggers
%         cfg.trialdef.eventvalue = 10; % Select trigger that is repeated only once in each trial
%         cfg.trialdef.prestim    = t_pre_stim;
%         cfg.trialdef.poststim   = t_post_stim;
%         cfg.event = evCorrected; % Event structure for the trials
%         cfg = ft_definetrial(cfg); % If your triggers are more complex, you might need to write your own functions
%         trl = cfg.trl; % all trials
%     
%         %% define LFP trials
%     
%         cfg.detrend = 'yes';
%         cfg.continuous = 'yes';
%         cfg.montage.tra      = bipolmat{v};
%         cfg.montage.labelnew = bipolabel{v};
%         cfg.montage.labelorg = labelorg{v};
%     
%         cfg.trl = trl
%         data = ft_preprocessing(cfg);
%     
%         onsets_for_spikes = trl(:,1:2)/hdr.Fs;
%     
%     
%         cfg.trl = trl(enc_onsets.eCorrRem,:);
%         onsets_for_eRspikes = cfg.trl(:,1:2)/hdr.Fs;
%         data_eR = ft_preprocessing(cfg);
%     
%         ekf = sort([enc_onsets.eMissed; enc_onsets.eCorrFam])
%     
%         cfg.trl = trl(ekf,:);
%         onsets_for_eKFspikes = cfg.trl(:,1:2)/hdr.Fs;
%         data_eKF = ft_preprocessing(cfg);
    
    
 %   file = ['pz', num2str(pz(isub)),'_uPreprodataSPC_prestim.mat'];
 %  save (file, 'hdr','evR','evCorrected','evClean','data','onsets_for_spikes', 'data_eR', 'data_eKF', 'onsets_for_eRspikes','onsets_for_eKFspikes');
%%

    load(fullfile(oripath,'data2github','SingleUnit',subjects{isub},['pz', num2str(pz(isub)),'_uPreprodataSPC_prestim.mat']))
       
%    downsample to 2kHz (following Jacobs et al., 2007)
    cfg             = [];
    cfg.resamplefs  = 2000;
    LFP_res         = ft_resampledata(cfg, data);
    LFP_res_eR      = ft_resampledata(cfg, data_eR);
    LFP_res_eKF      = ft_resampledata(cfg, data_eKF);
    
    
%    filter raw data for plotting
    cfg             = [];
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 40;
    LFP_4plot       = ft_preprocessing(cfg, LFP_res);
    LFP_4plot_eR       = ft_preprocessing(cfg, LFP_res_eR);
    LFP_4plot_eKF       = ft_preprocessing(cfg, LFP_res_eKF);
    
    freq_band = [3,12];
    fb = 1;
    for trl = 1:size(LFP_res.trial,2)
        signal{trl}          = LFP_res.trial{trl} - nanmean(LFP_res.trial{trl}); % center signal so that hilbert works properly
        filtfreqbounds  = [0, (1 - trans_width) * freq_band(fb,1), freq_band(fb,1), freq_band(fb,2), freq_band(fb,2) * (1 + trans_width), LFP_res.fsample/2] / (LFP_res.fsample / 2);
        filt_order      = round(2 * (LFP_res.fsample / freq_band(fb,1)));
        filterweights   = firls(filt_order, filtfreqbounds, idealresponse);
        filtered_signal{trl} = filtfilt(filterweights, 1, signal{trl}(1,:)); %when there are two bipolar channels I'm selecting the most lateral see channels definition
    end
      
%    bandpass filter raw data
    cfg             = [];
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [freq_band(fb,1) freq_band(fb,2)]; %here we bandpss filter between 3 and 20 Hz
    cfg.bpfiltord = 3; 
    LFP_bp       = ft_preprocessing(cfg, LFP_res);
    LFP_bp_eR    = ft_preprocessing(cfg, LFP_res_eR);
    LFP_bp_eKF   = ft_preprocessing(cfg, LFP_res_eKF);
    
%    hilbert the filtered signal
    for trl = 1:size(filtered_signal,2)
        temphilbert{trl}     = hilbert(filtered_signal{1,trl});
        anglehilbert{trl}    = angle(temphilbert{1,trl});
        bandphases{trl}      = anglehilbert{1,trl};
    end
        
%% wires to investigate
    
    %get available microwires
    wires   = dir(fullfile(spike_path, subjects{isub}, 'u*'));
    

    %% loop through wires
    for iwire = 1:size(wires, 1);
        
        fprintf('\n\nWire: %s.\n', wires(iwire).name);
        
        % load spiking data
        
        %load wave-clus output (for spike-depiction etc.)
        try
            t   = load(fullfile(wires(iwire).folder, wires(iwire).name, 'times_data.mat'));
        catch
            fprintf('No wave-clus for this wire.\n');
            continue;
        end
        
        
        
        % process spikes and relate them to the LFP        
        
        % loop through clusters
        for iclus = 1;%change this value to 2 to explore cluster 2:max(t.cluster_class(:, 1)) % do not analyze "0"-clusters at all (determined as "rest" by waveclus)
            
            fprintf('\t\tCluster: %d.\n', iclus);
            
            
            %get only data for this cluster (cave: you are analyzing more data points than when you are incorporating behavioral time)
            thiscluster{iwire}     = t.cluster_class(t.cluster_class(:, 1) == iclus, :); % cluster-number, microtime, behavioral time
            thisspike{iwire}       = t.spikes(t.cluster_class(:, 1) == iclus, :);
            
            thisclustersec{iwire} = thiscluster{iwire}./1000; %only active neurons with time in sec
            clustercllassec{iwire} = t.cluster_class./1000 % all spikes 1 and 0 with time in sec
            
            
            % associate spikes with theta phases
            
            onsets_for_spikes_align = (onsets_for_spikes(:,:) - onsets_for_spikes(1,1))
            spiketimeExp{iwire} = (thisclustersec{iwire}(:,2)- onsets_for_spikes(1,1));% align to the behavioural onset
            rawtodelete = spiketimeExp{iwire} <= 0; % select only spikes that occur during the behavioural exp
            spiketimeExp{iwire}(rawtodelete) = []
            thisspike{iwire}(rawtodelete,:) = []  
            
%% convert in ms
            spiketimeExp_ms{iwire} = spiketimeExp{iwire}.*1000;
            
           %% ISI refractoryness
           % = percentage of ISIs < 3ms
           clear ISI
           clear nISI
           ISI                 = diff(spiketimeExp_ms{iwire}(:, 1)); % inter-spike-intervals (ms)
           nISI                = size(ISI, 1); % number of ISIs
           percISIlessThan3ms{iwire}  = 100 * sum(ISI < 3) / nISI; % (%)
           %%
           [rows,cols,vals]= find(ISI < 3)
           isiTooShort{iwire} = rows
           spiketimeExp{iwire}(isiTooShort{iwire}) = []
           isi{iwire} = ISI;
           nisi{iwire}= nISI;
           
            %correct from the invalid trials due to trigger repetition
            
            load(fullfile(oripath,'data2github','SingleUnit',subjects{isub},['pz', num2str(pz(isub)),'m_data.mat']))
            tmp= indInvalidTrialsTimestamps;
            load(fullfile(oripath,'data2github','SingleUnit',subjects{isub},['pz', num2str(pz(isub)),'u_data.mat']))
            invtrl = unique([tmp; indInvalidTrialsTimestamps]);
            
            idp = LFP_res.time{1}(1,2) - LFP_res.time{1}(1,1);
            
            for trl = 2:size(LFP_res.trial,2)
            timetrlall(1,:) = LFP_res.time{1}(1,:);
            timetrlall(trl,:) = onsets_for_spikes_align(trl,1):idp:onsets_for_spikes_align(trl,2);
            end
            
            for trl=1:size(LFP_res.trial,2)
                bSpike{iwire}(trl,:)   	= hist(spiketimeExp{iwire},timetrlall(trl,:)) > 0;
                fprintf('Number of spikes: %d. Number of LFP time bins with >= 1 spike: %d. Difference = %d.\n', size(spiketimeExp{iwire}, 1), sum(bSpike{iwire}(trl,:)), size(spiketimeExp{iwire}, 1) - sum(bSpike{iwire}(trl,:)));
                [i,j,s] = find(bSpike{iwire}(trl,:) == 1);
                SpkLFP{iwire}(trl,:) = sum(s); % total nb of spike per trial
            end
 
 %           correct for spikes at 1 and 8000
            for trl=1:size(LFP_res.trial,2)
                if bSpike{iwire}(trl,1) == 1
                    bSpike{iwire}(trl,1) = 0
                else
                    bSpike{iwire}(trl,1) = bSpike{iwire}(trl,1)
                end
                if bSpike{iwire}(trl,end) == 1
                    bSpike{iwire}(trl,end) = 0
                else
                    bSpike{iwire}(trl,end) = bSpike{iwire}(trl,end)
                end
                fprintf('Number of spikes: %d. Number of LFP time bins with >= 1 spike: %d. Difference = %d.\n', size(spiketimeExp{iwire}, 1), sum(bSpike{iwire}(trl,:)), size(spiketimeExp{iwire}, 1) - sum(bSpike{iwire}(trl,:)));
                [i,j,s] = find(bSpike{iwire}(trl,:) == 1);
                SpkLFP_corr{iwire}(trl,:) = sum(s);
            end
            
%            get theta phases associated with the spike-time-bins
            
            for trl= 1:size(LFP_res.trial,2)
                sqallt = squeeze(bSpike{iwire}(trl,:));
                bpallt= squeeze(bandphases(1,trl));
                bphasesallt = bpallt{1,:};
                thetaPhasesallt{iwire}{trl} = bphasesallt(sqallt);
                binnedThetaPhasesallt{iwire}{trl}   = histcounts(thetaPhasesallt{iwire}{trl}, angle_bins);
                sqallt= []
                bpallt = []
                bphasesallt = []
            end
            
                 
            % select condition eR
            er = enc_onsets.eCorrRem
            % select condition eKF
            ekf = sort([enc_onsets.eMissed; enc_onsets.eCorrFam])
            
            er2rj = find(ismember(er,invtrl));
            er(er2rj) = []
                        
            ekf2rj = find(ismember(ekf,invtrl));
            ekf(ekf2rj) = []
            %% select data per each condition
            
            bSpikeeR{iwire} = bSpike{iwire}(er,:);
            bSpikeeKF{iwire} = bSpike{iwire}(ekf,:);
            
           
        end
    end
    
    cd (['C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\Prestim\' bins '\' cluster '\' direction])
    name = ['freq',num2str(freq_band(fb,1)),'to',num2str(freq_band(fb,2))]
    if ~exist(name)
        mkdir(name)
    end
    cd(name)
    
    file = ['Z',subjects{isub},'_data.mat'];
    save(file,'bSpikeeR', 'bSpikeeKF','onsets_for_spikes_align','spiketimeExp'); 


    close all
    clear binnedThetaPhaseseR
    clear binnedThetaPhaseseKF
    clear filtered_signal_eR
    clear filtered_signal_eKF
    clear bandphases_eR
    clear bandphases_eKF
    clear anglehilbert_eR
    clear anglehilbert_eKF
    clear AlleR
    clear AlleKF
    clear bSpikeeR
    clear bSpikeeKF
    clear trl_results
    clear trl_z_surro
    clear SPCeR
    clear SPCeKF
end


