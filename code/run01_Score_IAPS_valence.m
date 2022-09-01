clear all
close all
clc

% This script reproduce analysis and Supplementary Fig2b (valence)

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
load(fullfile(oripath,'data2github','Behaviour','all.txt')) %ALL STIMULI 1st 80 emotional; IAPS codes
load(fullfile(oripath,'data2github','Behaviour','arousal.mat'));% iaps number all, valence mean, valence sd, arousalmn, arousal sd. 

pict = (arousal(:,1))
vmn  = (arousal(:,2))
vsd  = (arousal(:,3))
amn  = (arousal(:,4))
asd  = (arousal(:,5))

[all_nw,ia,ib] = intersect(all, pict,'rows', 'stable')


IAPS_valence = arousal(ib,2)
IAPS_valence_sd = arousal(ib,3)

IAPS_select_valence = [all_nw, IAPS_valence,IAPS_valence_sd]
%%

subjects=[2 4 6 13 15 16 21 25 27 32 33 34 1 2 5 6 8 10]
mlist={'s2', 's4', 's6', 's13', 's15', 's16', 's21', 's25', 's27','s32','s33','s34','sz1','sz2','sz5','sz6','sz8','sz10'};

%%
numberofsubjects = size(subjects,2);
RT_IAPS_Enc = zeros(length(subjects),4);

index=0;
for sub=subjects
    index=index+1;
        if index<13;% these are the Ruber patients
        
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
    
    b= 500
    exclude = what_encoded < b 
    what_encoded(exclude) = []
   
   [what_encoded_e,ia_e,ib_e] = intersect(what_encoded,IAPS_select_valence(1:80));
   [what_encoded_n,ia_n,ib_n] = intersect(what_encoded,IAPS_select_valence(81:end));% 
   
    IAPS_valence_E = IAPS_select_valence(ib_e,2)
    IAPS_valence_slct_E = [what_encoded_e, IAPS_valence_E]
  
    IAPS_valence_N = IAPS_select_valence(80 + ib_n,2)
    IAPS_valence_slct_N = [what_encoded_n, IAPS_valence_N]
    
    e_pics_enc = find(ismember(what_encoded,what_encoded_e));%gives position of emotional pics in encoding presentation (40x1)
    n_pics_enc = find(ismember(what_encoded,what_encoded_n));%gives position of neutral pics in encoding presentation (80x1)
    
    e_pics_enc_ar = [e_pics_enc,IAPS_valence_E]
    n_pics_enc_ar = [n_pics_enc,IAPS_valence_N]
    
    e_items = find(IAPS_enc<41); %= e_words_enc but as a vector 1x40
    corr_e_items = setdiff(e_items,no_response_enc);
    n_items = find(IAPS_enc>40); %= n_words_enc but as a vector 1x80
    corr_n_items = setdiff(n_items,no_response_enc);

    enc_onsets = struct('e_items_corr',corr_e_items,'n_items_corr',corr_n_items); %e_ and n_items without not responded ones 


    %%%%RECOGNITION%%%%

       
    cd ../Rec
    load log.res
    response = log(:,1);

    load IAPS_rec
    what_recognised = all(IAPS_rec);

    b= 500
    exclude = what_recognised < b 
    what_recognised(exclude) = []
    
    
    old_E = find(ismember(what_recognised,what_encoded_e)); % all the old emotional pictures, listed by number of position they were shown
    what_recognised_e = find(ismember(what_recognised,all_nw(1:80))); % all emotional items shown (80), listed by number they were shown at
    new_E = setdiff(what_recognised_e,old_E); % all new emotional items (40) not shown the previous day, listed by number they were shown at 

    old_N = find(ismember(what_recognised,what_encoded_n));
    what_recognised_n = find(ismember(what_recognised,all_nw(81:end)));
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
    
    [eR,ia_re, ib_re] = intersect(old_eCorrRemID,IAPS_valence_slct_E)
    eR_ar = IAPS_valence_slct_E (ib_re,2)
    eR_valence = [eR eR_ar]
    
    old_eCorrFamID     = what_recognised(rec_onsets.eCorrFam);
    [eK,ia_ke, ib_ke] = intersect(old_eCorrFamID,IAPS_valence_slct_E)
    eK_ar = IAPS_valence_slct_E (ib_ke,2)
    eK_valence = [eK eK_ar]
    
    old_eMissedID  = what_recognised(rec_onsets.eMissed);
    [eMiss,ia_misse, ib_misse] = intersect(old_eMissedID,IAPS_valence_slct_E)
    eMiss_ar = IAPS_valence_slct_E (ib_misse,2)
    eMIss_valence = [eMiss eMiss_ar]
    
    old_eNotRespID = what_recognised(rec_onsets.eNotResp);
    [eNotResp,ia_eNotResp, ib_eNotResp] = intersect(old_eNotRespID,IAPS_valence_slct_E)
    eNotResp_ar = IAPS_valence_slct_E (ib_eNotResp,2)
    eNotResp_valence = [eNotResp eNotResp_ar]

     old_nCorrRemID     = what_recognised(rec_onsets.nCorrRem); %gets the stimulus number, IAPS reference number of the picture
    
    [nR,ia_rn, ib_rn] = intersect(old_nCorrRemID,IAPS_valence_slct_N)
    nR_ar = IAPS_valence_slct_N (ib_rn,2)
    nR_valence = [nR nR_ar]
    
    old_nCorrFamID     = what_recognised(rec_onsets.nCorrFam);
    [nK,ia_kn, ib_kn] = intersect(old_nCorrFamID,IAPS_valence_slct_N)
    nK_ar = IAPS_valence_slct_N (ib_kn,2)
    nK_valence = [nK nK_ar]
    
    old_nMissedID  = what_recognised(rec_onsets.nMissed);
    [nMiss,ia_missn, ib_missn] = intersect(old_nMissedID,IAPS_valence_slct_N)
    nMiss_ar = IAPS_valence_slct_N (ib_missn,2)
    nMIss_valence = [nMiss nMiss_ar]
    
    old_nNotRespID = what_recognised(rec_onsets.nNotResp);
    [nNotResp,ia_nNotResp, ib_nNotResp] = intersect(old_nNotRespID,IAPS_valence_slct_N)
    nNotResp_ar = IAPS_valence_slct_N (ib_nNotResp,2)
    nNotResp_valence = [nNotResp nNotResp_ar]

% plot
mvalence = [mean(eR_valence(:,2)), mean(nR_valence(:,2)), mean(eK_valence(:,2)), mean(nK_valence(:,2)), mean(eMIss_valence(:,2)),  mean(nMIss_valence(:,2)), mean(eNotResp_valence(:,2)), mean(nNotResp_valence(:,2))]

figure
bar(mvalence)
hold
title('valence')
set(gca,'XTick',[1:8])
set(gca,'XTickLabel',{'eR'; 'nR'; 'eK'; 'nK'; 'eMiss'; 'nMiss'; 'eNotresp'; 'nNotResp'})
ylabel('valence')

group_mvalence(index,:) = mvalence

end


%% load the data and plot the memory performance as function of valence
clear all
close all
clc

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils','beeswarm-master'));
load(fullfile(oripath,'data2github','Behaviour','group_valence_c1_c2.mat'));

mat = group_mvalence(:,[1:6])
figure('Position',[ 294   313   734   619],'Color','w');
ster = sterror(mat);
h = barwitherr(ster,nanmean(mat));
set(h(1),'FaceColor',[0.5 0.5 0.5]);
set(h(1),'EdgeColor',[1 1 1]);
set(gca,'XTickLabel',{'eR'; 'nR'; 'eK'; 'nK'; 'eF'; 'nF'}) % 'eNotresp'; 'nNotResp'
ylabel('valence');
hold;
plot(1.2,mat(:,1),'ko','MarkerFaceColor','k');
plot(2.2,mat(:,2),'ko','MarkerFaceColor','k');
plot(3.2,mat(:,3),'ko','MarkerFaceColor','k');
plot(4.2,mat(:,4),'ko','MarkerFaceColor','k');
plot(5.2,mat(:,5),'ko','MarkerFaceColor','k');
plot(6.2,mat(:,6),'ko','MarkerFaceColor','k');

group_mvalence = group_mvalence(:,[1:6]);

%save('group_valence_c1_c2.mat', 'group_mvalence');
