%%% Intracranial Patients
%%% IAPS Memory Test
%%% Day 1 Encode, Day 2 Recognition
%%% In recog test, R 97 - Left
%%%                K 100 - Down
%%%                N 98 - Right
%%% For encoding, which_e/n determines which of the pool of stimuli will be
%%% used for this part. which_e/n is used to write "encode_list"
%%% IAPS_enc is a randomisation of these stimuli (cogent uses IAPS_enc to
%%% index encode_list during encoding)

clear all
close all
clc

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','Behaviour','all.txt')) %ALL STIMULI 1st 80 emotional; IAPS codes
 
subjects=[2 4 6 13 15 16 21 25 27 32 33 34 1 2 5 6 8 10]
mlist={'s2', 's4', 's6', 's13', 's15', 's16', 's21', 's25', 's27','s32','s33','s34','sz1','sz2','sz5','sz6','sz8','sz10'};

numberofsubjects = size(subjects,2);
RT_IAPS_Enc = zeros(numberofsubjects,4);

index=0;
for sub=subjects
    index=index+1;
        if index<13;% these are Ruber patients
        
        %%%%ENCODING%%%%
        sdir = ['~\DataRuber\Patient',num2str(sub),'\Enc'];
        cd(sdir)
        
        else % load Zurich patients
        sdir=['~\DataZurich\Patient',num2str(sub),'\Enc'];
        cd(sdir)
    end
    %pics that subject did not respond to
    load log.res  %will need to exclude stimuli that patients do not respond to
    response_enc = log(:,1); no_response_enc = find(response_enc==0);
        
    load encode_list.txt %list of 120 IAPS pic numbers
    load IAPS_enc %% IAPS_enc is length 120; E items are 1:40, N from 41:120... which pertains to their position in Encode_List.txt
    
    
    what_encoded = encode_list(IAPS_enc); %GIVES THE EXACT ORDER IMAGES DISPLAYED
    
    what_encoded_e = intersect(what_encoded,all(1:80));
    what_encoded_n = intersect(what_encoded,all(81:end));
    
    e_pics_enc = find(ismember(what_encoded,what_encoded_e));%gives position of emotional pics in encoding presentation (40x1)
    n_pics_enc = find(ismember(what_encoded,what_encoded_n));%gives position of neutral pics in encoding presentation (80x1)
    
    e_items = find(IAPS_enc<41); %= e_words_enc but as a vector 1x40
    corr_e_items = setdiff(e_items,no_response_enc);
    n_items = find(IAPS_enc>40); %= n_words_enc but as a vector 1x80
    corr_n_items = setdiff(n_items,no_response_enc);
    
    
    %%%%RECOGNITION%%%%
        
    cd ../Rec
    
    load log.res
    response = log(:,1);
    
    
    load IAPS_rec
    what_recognised = all(IAPS_rec);
    
    old_E = find(ismember(what_recognised,what_encoded_e)); % all the old emotional pictures, listed by number of position they were shown
    what_recognised_e = find(ismember(what_recognised,all(1:80))); % all emotional items shown (80), listed by number they were shown at
    new_E = setdiff(what_recognised_e,old_E); % all new emotional items (40) not shown the previous day, listed by number they were shown at
    
    old_N = find(ismember(what_recognised,what_encoded_n));
    what_recognised_n = find(ismember(what_recognised,all(81:end)));
    new_N = setdiff(what_recognised_n,old_N);
    
    rec_onsets = struct('Old_e_items',old_E,'Old_n_items',old_N,'New_e_items',new_E,'New_n_items',new_N);
    rec_onsets.eCorrRem    =    old_E(find(response(old_E)==97));% correctly remembered
    rec_onsets.eCorrFam    =    old_E(find(response(old_E)==100));% correctly familiar with (known)
    rec_onsets.eMissed =    old_E(find(response(old_E)==98));% already seen but not recognised
    rec_onsets.eCorrNotRem   =    new_E(find(response(new_E)==98)); % never seen and not recognised
    rec_onsets.eFalseRem =    new_E(find(response(new_E)==97)); % wrongly remembered
    rec_onsets.eFalseFam =    new_E(find(response(new_E)==100));% wrongly familiar with
    rec_onsets.eNotResp = vertcat(old_E(find(response(old_E)==0)), new_E(find(response(new_E)==0)));  % No response in Recognition Task
    
    rec_onsets.nCorrRem    =    old_N(find(response(old_N)==97));
    rec_onsets.nCorrFam    =    old_N(find(response(old_N)==100));
    rec_onsets.nMissed =    old_N(find(response(old_N)==98));
    rec_onsets.nCorrNotRem   =    new_N(find(response(new_N)==98));
    rec_onsets.nFalseRem =    new_N(find(response(new_N)==97));
    rec_onsets.nFalseFam =    new_N(find(response(new_N)==100));
    rec_onsets.nNotResp = vertcat(old_N(find(response(old_N)==0)), new_N(find(response(new_N)==0)));
    
    
    old_eCorrRemID     = what_recognised(rec_onsets.eCorrRem); %gets the stimulus number, IAPS reference number of the picture
    old_eCorrFamID     = what_recognised(rec_onsets.eCorrFam);
    old_eMissedID  = what_recognised(rec_onsets.eMissed);
    old_eNotRespID = what_recognised(rec_onsets.eNotResp);
    
    old_nCorrRemID     = what_recognised(rec_onsets.nCorrRem);
    old_nCorrFamID     = what_recognised(rec_onsets.nCorrFam);
    old_nMissedID  = what_recognised(rec_onsets.nMissed);
    old_nNotRespID = what_recognised(rec_onsets.nNotResp);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%Now go back to Encoding to do subsequent memory analysis%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    enc_onsets = struct('e_items_all',e_pics_enc,'n_items_all',n_pics_enc); 
    
    enc_onsets.eCorrRem = find(ismember(what_encoded,old_eCorrRemID));
    enc_onsets.eCorrFam = find(ismember(what_encoded,old_eCorrFamID));
    enc_onsets.eMissed = find(ismember(what_encoded,old_eMissedID));
    enc_onsets.eNotResp = find(ismember(what_encoded,old_eNotRespID));
    
    enc_onsets.nCorrRem = find(ismember(what_encoded,old_nCorrRemID));
    enc_onsets.nCorrFam = find(ismember(what_encoded,old_nCorrFamID));
    enc_onsets.nMissed = find(ismember(what_encoded,old_nMissedID));
    enc_onsets.nNotResp = find(ismember(what_encoded,old_nNotRespID));
    
    %% plots
    % Now do memory performance analysis
    
    recog(index,:) = [length(find(response(old_E)==97))/length(old_E); %E R Hits
        length(find(response(new_E)==97))/length(new_E); %E R FAs
        length(find(response(old_E)==100))/length(old_E); %E K Hits
        length(find(response(new_E)==100))/length(new_E); %E K FAs
        length(find(response(old_N)==97))/length(old_N); %N R Hits
        length(find(response(new_N)==97))/length(new_N); %N R FAs
        length(find(response(old_N)==100))/length(old_N); %N K Hits
        length(find(response(new_N)==100))/length(new_N); %N K FAs
        length(find(response(old_E)==98))/length(old_E); %E Miss
        length(find(response(old_N)==98))/length(old_N); %N Miss
        length(find(response(new_E)==98))/length(new_E); %E Crj
        length(find(response(new_N)==98))/length(new_N);] %N Crj
    
end  

%% load the data and plot the group behavioural results
clear all
close all
clc

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'));
load(fullfile(oripath,'data2github','Behaviour','beh_cohort1.mat'));
load(fullfile(oripath,'data2github','Behaviour','beh_cohort2.mat'));

mat= [cohort1; cohort2]

%% Reproduce Figure 1b
x = [ones(18,1) ones(18,1)*2 ones(18,1)*3 ones(18,1)*4 ones(18,1)*5 ones(18,1)*6 ones(18,1)*7 ones(18,1)*8];
y = [mat(:,1) mat(:,2) mat(:,3) mat(:,4) mat(:,5) mat(:,6) mat(:,7) mat(:,8)];

figure;
beeswarm(x(:),y(:),'sort_style','up','dot_size',4,'overlay_style','sd','colormap',[1 0 0; 0.7 0 0; 0 0 0; 0.2 0.2 0.2; 0 0 1; 0 0 0.7;0.5 0.5 0.5; 0.7 0.7 0.7])
ylim([-9 100]);
xticklabels({'eR','eRFA','eK','eKFA','nR','nRFA','nK','nKFA'})
%% stat 
load(fullfile(oripath,'data2github','Behaviour','behav_mat_spss_c1c2.mat'));% data from n=18 patients from cohort1 n=13 and cohort 2 n=5, response types: Re Ke Rn Kn RFase RFAsn

%% calculate mean, sem
clear mat
mat = behav_cohort1cohort2_spss

dfase = mean(mat(:,5))
dfasn = mean(mat(:,6))
semfase = sem(mat(:,5))
semfasn = sem(mat(:,6))

dke = mean(mat(:,2))
dkn = mean(mat(:,4))
semke = sem(mat(:,2))
semkn = sem(mat(:,4))

appmat_fas = [mat(:,6); mat(:,5)];
stdfas = std(appmat_fas)

appmat_k = [mat(:,4); mat(:,2)];
stdk = std(appmat_k)

dfas = (dfase-dfasn)/stdfas
dk= (dke-dkn)/stdk