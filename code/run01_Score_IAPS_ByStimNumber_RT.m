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

%% This script produce the results in Group18_IAPS_RTs_c1c2.mat and Group18_IAPS_RTs_rec_c1c2.mat which are reported in Supplementary Table 3 and 4

clear all
close all
clc

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','Behaviour','all.txt')) %ALL STIMULI 1st 80 emotional; IAPS codes

subjects=[2 4 6 13 15 16 21 25 27 32 33 34 1 2 5 6 8 10];
mlist={'s2', 's4', 's6', 's13', 's15', 's16', 's21', 's25','s27','s32','s33','s34','sz1','sz2','sz5','sz6','sz8','sz10'};

numberofsubjects = size(subjects,2); 
RT_IAPS_Enc = zeros(numberofsubjects,4); 

index=0;
for sub=subjects
    index=index+1;
        if index<13
        
        %%%%ENCODING%%%%
        sdir = ['~\DataRuber\Patient',num2str(sub),'\Enc'];
        cd (sdir)
        
    else
        sdir=['~\DataZurich\Patient',num2str(sub),'\Enc'];
        cd(sdir)
    end

    %%%%ENCODING%%%%

   
    %pics that subject did not respond to
    load log.res  %will need to exclude stimuli that patients do not respond to
    response_enc = log(:,1); no_response_enc = find(response_enc==0); 
    RT_enc = log(:,2)
    
    load IAPS_enc %% IAPS_enc is length 120; E items are 1:40, N from 41:120... which pertains to their position in Encode_List.txt
    emoRT = find(IAPS_enc<41); 
    emoRT = RT_enc(emoRT);
    emoRT = emoRT(find(emoRT));
    neuRT = find(IAPS_enc>40); 
    neuRT = RT_enc(neuRT);
    neuRT = neuRT(find(neuRT));
    emosem = std(emoRT)/sqrt(length(emoRT));
    neusem = std(neuRT)/sqrt(length(neuRT));
   
    
    RT_IAPS_Enc(index,:) = round([mean(emoRT) mean(neuRT) emosem neusem ]);
    
end

meanRT = mean(RT_IAPS_Enc(:,[1:2]));
semRT = std(RT_IAPS_Enc(:,[1:2]))/sqrt(length(subjects))

% save('Group18_IAPS_RTs_c1c2.mat','RT_IAPS_Enc')
%% 
clear all
close all
clc

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','Behaviour','all.txt')) %ALL STIMULI 1st 80 emotional; IAPS codes

subjects=[2 4 6 13 15 16 21 25 27 32 33 34 1 2 5 6 8 10];
mlist={'s2', 's4', 's6', 's13', 's15', 's16', 's21', 's25','s27','s32','s33','s34','sz1','sz2','sz5','sz6','sz8','sz10'};


numberofsubjects = size(subjects,2); 
RT_IAPS_Enc = zeros(numberofsubjects,4); 

index=0;
for sub=subjects
    index=index+1;
    if index<13
        
        %%%%ENCODING%%%%
        sdir = ['~\DataRuber\Patient',num2str(sub),'\Enc'];
        cd (sdir)
        
    else
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
    
    cd ../Rec
     
    %pics that subject did not respond to
    load log.res
    response = log(:,1);
    RT_rec = log(:,2)
   
    load IAPS_rec
    what_recognised = all(IAPS_rec);

    old_E = find(ismember(what_recognised,what_encoded_e)); % all the old emotional pictures, listed by number of position they were shown
    what_recognised_e = find(ismember(what_recognised,all(1:80))); % all emotional items shown (80), listed by number they were shown at
    new_E = setdiff(what_recognised_e,old_E); % all new emotional items (40) not shown the previous day, listed by number they were shown at 

    old_N = find(ismember(what_recognised,what_encoded_n));
    what_recognised_n = find(ismember(what_recognised,all(81:end)));
    new_N = setdiff(what_recognised_n,old_N);
    
    %%
    eCorrRemRT =   old_E(find(response(old_E)==97));
    eCorrRemRT = RT_rec(eCorrRemRT);
    
    eCorrFamRT =   old_E(find(response(old_E)==100));
    eCorrFamRT = RT_rec(eCorrFamRT);
    
    eMissedRT = old_E(find(response(old_E)==98));
    eMissedRT = RT_rec(eMissedRT);
    
    eCorrNotRemRT = new_E(find(response(new_E)==98));
    eCorrNotRemRT = RT_rec(eCorrNotRemRT);
    
    eFalseRemRT = new_E(find(response(new_E)==97)); 
    eFalseRemRT = RT_rec(eFalseRemRT);
    
    eFalseFamRT = new_E(find(response(new_E)==100));
    eFalseFamRT = RT_rec(eFalseFamRT);
    %%
    nCorrRemRT =   old_N(find(response(old_N)==97));
    nCorrRemRT = RT_rec(nCorrRemRT);
    
    nCorrFamRT =   old_N(find(response(old_N)==100));
    nCorrFamRT = RT_rec(nCorrFamRT);
    
    nMissedRT =   old_N(find(response(old_N)==98));
    nMissedRT = RT_rec(nMissedRT);
    
    nCorrNotRemRT = new_N(find(response(new_N)==98));
    nCorrNotRemRT = RT_rec(nCorrNotRemRT);
    
    nFalseRemRT = new_N(find(response(new_N)==97)); 
    nFalseRemRT = RT_rec(nFalseRemRT);
    
    nFalseFamRT = new_N(find(response(new_N)==100));
    nFalseFamRT = RT_rec(nFalseFamRT);
    
    eR = std(eCorrRemRT)/sqrt(length(eCorrRemRT));
    eFasR = std(eFalseRemRT)/sqrt(length(eFalseRemRT));
    eK = std(eCorrFamRT)/sqrt(length(eCorrFamRT));
    eFasK = std(eFalseFamRT)/sqrt(length(eFalseFamRT));
    eMiss = std(eMissedRT)/sqrt(length(eMissedRT));
    eCRj = std(eCorrNotRemRT)/sqrt(length(eCorrNotRemRT));
    
    nR = std(nCorrRemRT)/sqrt(length(nCorrRemRT));
    nFasR = std(nFalseRemRT)/sqrt(length(nFalseRemRT));
    nK = std(nCorrFamRT)/sqrt(length(nCorrFamRT));
    nFasK = std(nFalseFamRT)/sqrt(length(nFalseFamRT));
    nMiss = std(nMissedRT)/sqrt(length(nMissedRT));
    nCRj = std(nCorrNotRemRT)/sqrt(length(nCorrNotRemRT));
    
  

    RT_IAPS_Rec_emo(index,:) = round([mean(eCorrRemRT) mean(eFalseRemRT) mean(eCorrFamRT) mean(eFalseFamRT) mean(eMissedRT) mean(eCorrNotRemRT) eR eFasR eK eFasK eMiss eCRj]);
    RT_IAPS_Rec_neu(index,:) = round([mean(nCorrRemRT) mean(nFalseRemRT) mean(nCorrFamRT) mean(nFalseFamRT) mean(nMissedRT) mean(nCorrNotRemRT) nR nFasR nK nFasK nMiss nCRj]); 
    
end

meanRT_emo = nanmean(RT_IAPS_Rec_emo(:,[1:6]))
semRT_emo = nanstd(RT_IAPS_Rec_emo(:,[1:6]))/sqrt(length(subjects))

meanRT_neu = nanmean(RT_IAPS_Rec_neu(:,[1:6]))
semRT_neu = nanstd(RT_IAPS_Rec_neu(:,[1:6]))/sqrt(length(subjects))

%save('Group18_IAPS_RTs_rec_c1c2.mat','RT_IAPS_Rec_emo', 'RT_IAPS_Rec_neu');

