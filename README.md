# CostaLozanoetal

This folder contains the code used to analyze the data included in the paper

Behavioural data analysis (01):
- run01_Score_IAPS_ByStimNumber_Behaviour.m. Code to process and analyze memory performance for each patient and for all trials types. The code reproduce Figure 1b.
- run01_Score_IAPS_ByStimNumber_RT.m. Code to process and analyze reaction time at encoding as reported in Supplementary Table 4 and 5.
- run01_Score_IAPS_arousal.m. Code to  reproduce analysis and Supplementary Fig. 3a.
- run01_Score_IAPS_valence.m. Code to  reproduce analysis and Supplementary Fig. 3b.

Preprocessing and time frequency analysis and cleaning of the electrophisiological data (02).  

- run02_Prepro_bipolar_Amy.m.Code to preprocess the amygdala data for Ruber patients
- run02_Prepro_bipolar_Amy_pz.m.Code to preprocess the amygdala data for Zurich patients
- run02_Prepro_bipolar_Hippo.m.Code to preprocess the hippocampal data for Ruber patients
- run02_Prepro_bipolar_Hippo_pz.m.Code to preprocess the hippocampal data for Zurich patients
- run02_Prepro_bipolar_ECPRh.m.Code to preprocess the Entorhinal (EC) and Perirhinal (PRh) data for Ruber patients
- run02_Prepro_bipolar_ECPRh_pz.m.Code to preprocess the Entorhinal (EC) and Perirhinal (PRh) data for Zurich patients
- run02_TFcleaning_Amy.m.Code to performe time frequency analysis in the amygdala and trial by trial visual inspection for rejecting artifacts 
- run02_TFcleaning_Hc.m.Code to performe time frequency analysis in the hippocampus and trial by trial visual inspection for rejecting artifacts 
- run02_TFcleaning_ECPRh.m.Code to performe time frequency analysis in the Entorhinal (EC) and Perirhinal (PRh) and trial by trial visual inspection for rejecting artifacts 

Time frequency results and statistics reported in the paper (03)
- run03_TFstat_Amyall reproduce the time frequency results and statistics for patients with electrodes implanted in the amygdala. It also produces Figure 1d,e,f and Supplementary Fig.5.
- run03_TFstat_Amy reproduce the time frequency results and statistics for patients with electrodes implanted in the amygdala and that also have electrodes in the hippocampus. It produces Supplementary Figure 4
- run03_TFstat_Hippo reproduce the time frequency results and statistics for patients with electrodes implanted in the hippocampus. It also produces Figure 1g,h,i and Supplementary Fig.5 and 7
- run03_TFstat_ECPRh reproduce the time frequency results and statistics for patients with electrodes implanted in the Entorhinal (EC) and Perirhinal (PRh). It also produces Supplementary Fig.8b,c,d

Connectivity (04)
- run04_coherence.m.Code to run coherence between stimulus-induced amygdala and hippocampal activity. The code produces Supplementary Fig.9
- run04_GrangerCausality.m. Code to run time frequency granger causality in the direction amygdala to hippocampus and hippocampus to amygdala as reported in Fig.2 a,b,c and Supplementary Fig.10
- run04_PAC_AMYAMY.m Code to compute the phase amplitude coupling (PAC) within the amygdala,Supplementary Fig 12,14
- run04_PAC_AMYAMY_subsampled.m Code to compute the phase amplitude coupling (PAC) within the amygdala controlling for number of trials, Supplementary Fig 12,14
- run04_PAC_AMYHC.m Code to compute the phase amplitude coupling (PAC) between amygdala phase and hippocampus amplitude, Figure 2d and Supplementary Fig. 11,13
- run04_PAC_AMYHC_subsampled.m Code to compute the phase amplitude coupling (PAC) between amygdala phase and hippocampus amplitude controlling for number of trials, Supplementary Fig. 11,13
- run04_PAC_HCHC.m Code to compute the phase amplitude coupling (PAC) within the hippocampus,Supplementary Fig 12,14
- run04_PAC_HCHC_subsampled.m Code to compute the phase amplitude coupling (PAC) within the hippocampus,Supplementary Fig 12,14 controlling for number of trials
- run04_PAC_HCAMY.m Code to compute the phase amplitude coupling (PAC) between hippocampus phase and amygdala amplitude, Supplementary Fig 11,13
- run04_PAC_H.mCAMY_subsampled.m Code to compute the phase amplitude coupling (PAC) between hippocampus phase and amygdala amplitude controlling for number of trials, Supplementary Fig 11,13
- run04_PACOi This script computes the phase to amplitude coupling opposition index (PACOi)
 between amygdala phase and hippocampus amplitude contacts between eR and eKF
 trials
- run04_PTA.m This script computes the peak trigger averages (PTA) within and between regions
- run04_corr_envgpeakslag_PACOi.m This script computes the correlation between PACOi and the amplitude envelope cross-correlation method. Fig 5b, Supplementary Fig. 20


Spike sorting and basic about units(05)
- run05_ExtractMicroData.m This script extracts the microwire data from the original files for use with wave-clus.
- run05_GetSpikesUsingWaveClus.m This script performs spike-detection and -clustering using wave-clus.
- run05_Basicaboutunits.m. This script compute the quality assessment of neuronal recordings and spike sorting collecting informations across units found in patient 5,6,8,10 
- run05_Basicaboutunits_special.m. This script do the same as before but for patient 1 and 2 where we corrected for an error due to a system data acquisition error
- run05_Plots_basicunits.m. This script reproduce Supplementary Fig. 21 a,b,c,d

SPC analysis (06)
- run06_SPC_amyhc.m. Code to performe spike field coherence analysis between amygdala theta and hippocampal spikes for pz5, 6, 8, 10
- run06_SPC_amyhc_special.m Code to performe spike field coherence analysis between amygdala theta and hippocampal spikes for pz1, 2. For two patients the LFP the LFP was extracted from amygdala microelectrode due to a system data acquisition error
- run06_Singlesubject_SPC_results.m  Code to put together spike field coherence results for each unit.
- run06_SPCGroup_results.m Code to plot spike field coherence results as average spikes of all units. The code reproduce Fig 3 b and Supplementary Fig. 24
- run06_Singlesubject_Rasterplots.m plot the event related raster plots for each unit for eR and eKF condition from -0.5 to 1.5 sec as in Supplementary Fig.22.
- run06_Oneexampletrl_SPC Code to produce supplementary Fig. 23, an example trial showing the occurrence of hippocampal spikes from one unit in relation with the ipsilateral amygdala theta phase
