
%% concatenate all units found from positive and negative cluster (1|2)
% clear all
% close all
% clc
% 
% bins = '20_bins' 
% 
% 
% load('pos_c2.mat')
% load('pos_c1.mat')
% load('neg_c1.mat')
% 
% eMa_all = [eMa1; eMa2; eMa3];
% nMa_all = [nMa1; nMa2; nMa3];
% 
% clear eMa
% clear nMa
% 
% eMa = eMa_all
% nMa = nMa_all
% 
% % 31 units for both conditions, phase realigned
% sumeMa = sum(eMa_all,2)
% sumnMa = sum(nMa_all,2)

%% this script produce Figure 3f, g  and Supplementary Fig 24
clear all
close all
clc

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','SingleUnit','Results_allunits.mat'))

pathfigs = '/figs/';% subpath to save figures and results
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','pac'))
addpath(fullfile(oripath,'Costalozanoetal','code','pacoi'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','circstat-matlab'))

angle_bins      = deg2rad(-180:18:180);
angle_centers  	= transpose(movmean(angle_bins, 2, 'endpoints', 'discard'));
x_dbin =angle_centers'

eMc = mean(eMa.*exp(sqrt(-1)*x_dbin),2); % Ema is spikes occuring in each bin for each unit
nMc = mean(nMa.*exp(sqrt(-1)*x_dbin),2);

E_ph = angle(mean(eMc./abs(eMc))); %preferential phase for eR average over unit
N_ph = angle(mean(nMc./abs(nMc)));

E_ab = abs(mean(eMc./abs(eMc)));
N_ab = abs(mean(nMc./abs(nMc)));

rE = mean(zscore(eMa,[],2));
rN = mean(zscore(nMa,[],2));

rErose = mean(eMa);
rNrose = mean(nMa);

rErosenorm_all = eMa./sum(eMa,2);
rNrosenorm_all = nMa./sum(nMa,2);

rErosenorm = mean(rErosenorm_all);
rNrosenorm = mean(rNrosenorm_all);
%% average of all units n=31

x = linspace(-0.5,0.5,20);
w=2*pi*1;
W=[w-w/8];
[param1,y_est1,y_resid1,err_rms1] = sinefit(rErose,x,W);
[param2,y_est2,y_resid2,err_rms2] = sinefit(rNrose,x,W);

h1 = Figure(1,'size',[120   120],'fontSizeScreen',10);
subplot(221);bar(x_dbin,(rErose),'FaceColor','none','EdgeColor',[1 0 0],'LineWidth',1);
subplot(221);hold all;bar(x_dbin,(rNrose),'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1);
subplot(221);plot(E_ph,5,'vr','MarkerSize',8,'MarkerFaceColor','r')
subplot(221);plot(N_ph,5,'vk','MarkerSize',8,'MarkerFaceColor','k')
xlabel('AMY_\theta phase (rad)');
ylabel('Average hippocampus spikes count');

subplot(221);plot(x_dbin,y_est1,'Color',[1 0 0],'LineWidth',2);
subplot(221);hold all;plot(x_dbin,y_est2,'Color',[0 0 0],'LineWidth',2);
xlim([-pi-0.1 pi+0.1]);
set(gca, 'box','off')
ll=legend({'eR','eKF','eR angle','eKF angle'},'location','northwest');
set(gca, 'box','off')
set(ll, 'box','off')

%% group stat
[h1,p1,ks2stat1] = kstest2(rErose, rNrose)
fprintf('%5f\n',p1)

%% permutation using Watson-Williams test
po = angle(eMc./abs(eMc)) - angle(nMc./abs(nMc));
    ang_eR = angle(eMc./abs(eMc));
    ang_eKF = angle(nMc./abs(nMc));
    
nReps = 10000;
perm = zeros(nReps,1);
[pval, table]=circ_wwtest(ang_eR,ang_eKF);
f_obs = table{2,5};
p_obs = pval;
npat=size(ang_eR,1);

mat=cat(1,ang_eR,ang_eKF);
for i=1:nReps
  [null,index] = sort(rand(size(mat)));
  sel1 = index(1:npat);
  sel2 = index(npat+1:end);
  [dum,table] = circ_wwtest(mat(sel1,1),mat(sel2,1));
  perm(i)=table{2,5};
end
p = (sum(perm>f_obs)+1)/(nReps+1);

figure;hist(perm,50);
ylim = get(gca,'YLim');
hold all;plot(f_obs*[1,1],ylim,'r-','LineWidth',2);
xlabel('Watson-Williams F value')
ylabel('counts')
%% single cell test using circ_kuipertest, values reported on supplementary Fig.24

for u=1:31
    [pval1,k1,K1] = circ_kuipertest(eMa(u,:),nMa(u,:)) 
    p_all(u,1) = pval1;
end
%%

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','SingleUnit','wave_pc1.mat'))
pc1 = thisspike_sj
load(fullfile(oripath,'data2github','SingleUnit','wave_nc1.mat'))
nc1 = thisspike_sj
load(fullfile(oripath,'data2github','SingleUnit','wave_pc2.mat'))
pc2 = thisspike_sj

allw = [pc1 nc1 pc2] 

count=0
for i=1:21
   
    for u = 1:size(allw{i},2)
        count=count+1
        wv{count} = allw{1,i}{1,u}
    end
end

%save ('n31_wave.mat', 'wv');
%% plot data
bins = '20_bins' 
angle_bins      = deg2rad(-180:18:180);
angle_centers  	= transpose(movmean(angle_bins, 2, 'endpoints', 'discard'));

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils','circstat-matlab'))
load(fullfile(oripath,'data2github','SingleUnit','20_binsFig_data.mat'))
load(fullfile(oripath,'data2github','SingleUnit','n31_wave.mat'))

FruHzn31 = FruHz([1:13,15:32],:)

%%
AlleMa = eMa
AllnMa = nMa

%% Figure 3g
clear dt
dt.figEx = [2 5 12 8]
dt.boxEx = [4 1.35 4.25 4.25]
dt.visible = 'on'
dt.angLabels{1} = '0°'
dt.angLabels{2} = '90°'
dt.angLabels{3} = '180°'
dt.angLabels{4}= '-90°'
dt.bReverseY = 1
dt.angleCenters = angle_centers
dt.angularRes = 20
dt.FR = sum(AlleMa)'
dt.FR1 = sum(AllnMa)'
dt.t.par.sr = 32000;
dt.nspk = sum(sum(AlleMa))
dt.nspk1 = sum(sum(AllnMa))

h= plotFRDir_group(dt)
% set(h, 'PaperPositionMode', 'auto');
% print(h, strcat('groupres_erekfv1'), '-dtiff', '-r450');

% save data to recreate figure
save('groupres_erekfv1');
close all
%% plot each unit both conditions
% Supplementary Figure 24

addpath ('C:\Users\manuela\Desktop\CTB\00000_SubmittingIAPS\CostaLozanoetal\code\utils')

for u=1:31
clear dt
dt.figEx = [2 5 12 8]
dt.boxEx = [4 1.35 4.25 4.25]
dt.visible = 'on'
dt.angLabels{1} = '0°'
dt.angLabels{2} = '90°'
dt.angLabels{3} = '180°'
dt.angLabels{4}= '-90°'
dt.bReverseY = 1
dt.angleCenters = angle_centers
dt.angularRes = 20
dt.FR = (AlleMa(u,:))';
dt.FR1 = (AllnMa(u,:))';
dt.waveEx = [1 5.5 2 2]
dt.t.par.sr = 32000;
dt.thisSpike = wv{1,u}
dt.nspk = sum(AlleMa(u,:))
dt.nspk1 = sum(AllnMa(u,:))


h= plotFRDir_group(dt)
set(h, 'PaperPositionMode', 'auto');
%print(h, strcat('eReKF_UnitsFRplot', num2str(u)), '-dtiff', '-r450');

%save data to recreate figure
%save(strcat('uniteReKF', num2str(u)));
close all
end

%% calculate the delay from the phase opposition
delay = ((1000/7.5)*1.32)/(2*pi)