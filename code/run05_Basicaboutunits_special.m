clear all; close all; clc;

% add fieldtrip
restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','eeg_toolbox_v1_3_2'))% kahana eeg toolbox
addpath(fullfile(oripath,'Costalozanoetal','code','utils','MatlabImportExport_v6.0.0'));
warning('off', 'all');% suppress all warnings

% select manually depending on wether you want to look at the data 
% in the positive or negative spikes folder , either cluster 1 or cluster 2

direction='posSpikes';% select manually
%direction='negSpikes';

cluster= 'Cluster1';% select manually
%cluster= 'Cluster2'

bins = '20_bins'

% direct to the folder with wave clus results
spike_path      = ['C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\SpikeExtraction_20210210\' direction];


% settings
angle_bins      = deg2rad(-180:18:180); % if 20 bins are choosen
angle_centers  	= transpose(movmean(angle_bins, 2, 'endpoints', 'discard'));
t_pre_stim = 0; % How many seconds before the trigger to take - scenes were presented for 500 ms + 3500 of iti
t_post_stim = 4; % How many seconds after the trigger to take


trans_width     = 0.15;
idealresponse   = [0, 0, 1, 1, 0, 0];
numSurrogates   = 1000;
cutoffHighPower = 0;

% analysis options
qRemoveSpikes           = 'NoremSpikes'; % question whether to remove spikes from the LFP before processing it: 'remSpikes' or 'notRemSpikes'
data4Rayleigh           = 'raw4Rayleigh'; % whether to normalize data before Rayleigh test: 'raw4Rayleigh' or 'norm4Rayleigh'
surroCreationType       = 'fullyRandom'; % how to create the surrogate data: 'fullyRandom' or 'closeAround'

% subjects list/ a special script is done for pz1 and pz2
subjects    = {...
    'Patient1'; ...
    'Patient2'; ...
    'Patient200'; ...
    };

pz = [1 2 200];
sj = {'pz1u', 'pz2u', 'pz200u'}

allRes          = [];
i=0
v=0
%% loop through subjects
for isub = 1:length(subjects)
    v=v+1;
    i=i+1;
    
    fprintf('\n\nSUBJECT: %s.\n', subjects{isub});
   
    
    %% get onsets
    load(fullfile(oripath,'data2github','SingleUnit','From0to4sec',['pz', num2str(pz(isub)),'_PreprodataSPC.mat']))
    
    %% wires to investigate
    
    % get available microwires
    wires   = dir(fullfile(spike_path, subjects{isub}, 'u*'));

    %% loop through wires
    for iwire = 1:size(wires, 1);
        
        fprintf('\n\nWire: %s.\n', wires(iwire).name);
        %% load spiking data
        
        % load wave-clus output (for spike-depiction etc.). Time data
        % contain action potentials 1 or 2 (not use 0) and the time in ms
        % at which they occur.
        try
            t   = load(fullfile(wires(iwire).folder, wires(iwire).name, 'times_data.mat'));
        catch
            fprintf('No wave-clus for this wire.\n');
            continue;
        end
        
        
        try
            raw   = load(fullfile(wires(iwire).folder, wires(iwire).name, 'data.mat'));
        catch
            fprintf('No wave-clus for this wire.\n');
            continue;
        end
        

        %% focus only on the raw data during the experiment
        dataexp_min = (size(data)/raw.sr)/60
        
        load(fullfile(oripath,'data2github','SingleUnit','From0to4sec',[sj{i}, '_data.mat']))

        timestamp = ((evCorrected([evCorrected.value] == 10)));
        begexp = (timestamp(1,1).sample/raw.sr)/60;
        endexp = (timestamp(1,end).sample/raw.sr)/60;
        beh = Onset_times(end,1)/60000
        
        datacut = raw.data(timestamp(1,1).sample:timestamp(1,end).sample);
        
        par                 = [];
        par.sr              = raw.sr;
        par.detect_fmin     = t.par.detect_fmin;
        par.detect_fmax     = t.par.detect_fmax;
        par.detect_order    = t.par.detect_order;
        xf_detect           = spike_detection_filter(datacut, par);
        
        % calculate STD of the filtered data
        STDnoise            = std(xf_detect);
        %% process spikes and relate them to the LFP     
        %% loop through clusters
        for iclus = 1:max(t.cluster_class(:, 1)) % do not analyze "0"-clusters at all (determined as "rest" by waveclus)
            
            fprintf('\t\tCluster: %d.\n', iclus);
            
            % get only data for this cluster (cave: you are analyzing more data points than when you are incorporating behavioral time)
            thiscluster{iwire}     = t.cluster_class(t.cluster_class(:, 1) == iclus, :); % cluster-number, microtime, behavioral time
            thisspike{iwire}       = t.spikes(t.cluster_class(:, 1) == iclus, :);
            
            thisclustersec{iwire} = thiscluster{iwire}./1000; %only active neurons with time in sec
            clustercllassec{iwire} = t.cluster_class./1000 % all spikes 1 and 0 with time in sec
            
            
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
           percISIlessThan3ms  = 100 * sum(ISI < 3) / nISI; % (%)
           %%
           [rows,cols,vals]= find(ISI < 3)
           isiTooShort{iwire} = rows
           spiketimeExp{iwire}(isiTooShort{iwire}) = []
           isi{iwire} = ISI;
           nisi{iwire}= nISI;
           percISItooshort{iwire}= percISIlessThan3ms;
           %%
           meanFR  	= size(spiketimeExp_ms{iwire}(:, 1),1) / (range(spiketimeExp_ms{iwire}(:, 1)) / 1000); % (Hz)          
           %% waveform peak SNR (Faraut et al., 2018)
            % (= ratio between the peak amplitude of the mean waveform and
            % the STD of the noise)  
            param.peakIdx    	= 20; % peak index of the waveform
            
            tspk = thisspike{iwire};
            peakAmpl   	= abs(mean(thisspike{iwire}(:, param.peakIdx)));
            peakSNR  	= peakAmpl / STDnoise;   
            
            %% collect information across units
            
            % this unit's results
            unitRes                     = [];
            unitRes.idx                 = [isub, iwire, iclus];
            unitRes.percISIlessThan3ms  = percISIlessThan3ms;
            unitRes.meanFR              = meanFR;
            unitRes.peakSNR             = peakSNR;
            
            % collapse across units
            allRes  = cat(1, allRes, unitRes);
            
        end
    end
    
end