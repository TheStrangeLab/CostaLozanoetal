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

direction='posSpikes';% posSpikes | negSpikes 
cluster= 'Cluster1';% Cluster1 | Cluster2

bins = '20_bins' 


% direct to the folder with wave clus results
spike_path      = ['C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\SpikeExtraction_20210210\' direction];

% subjects list
subjects    = {...
    'Patient1'; ...
    'Patient2'; ...
    'Patient200'; ...
    'Patient5'; ...
    'Patient6'; ...
    'Patient8'; ...
    'Patient10'; ...
    };

angle_bins      = deg2rad(-180:18:180);%deg2rad(-180:3:180);
angle_centers  	= transpose(movmean(angle_bins, 2, 'endpoints', 'discard'));


v=0


iclus = 1% change this value when you want to look at cluster 2
for isub = 1:length(subjects)
    
    freq_band = [3,12];
    %for
    fb = 1
    v=v+1;
    cd (['C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\Prestim\' bins '\' cluster '\' direction]);
    name = ['freq',num2str(freq_band(fb,1)),'to',num2str(freq_band(fb,2))]
    cd (name)
    load(['Z',subjects{isub},'_data']);
    
    for iwire = 1:size(bSpikeeR,2)
       
        if size(bSpikeeR{1,iwire},1 ~= 0)
        for trl = 1:size(bSpikeeR{1,iwire},1)
        eR{isub}{iwire}(trl,:) = bSpikeeR{1,iwire}(trl,:)
        end
        else
        fprintf('No wave-clus for this wire.\n');
        end
    end
    
    for iwire = 1:size(bSpikeeKF,2)
       
        if size(bSpikeeKF{1,iwire},1 ~= 0)
        for trl = 1:size(bSpikeeKF{1,iwire},1)
        eKF{isub}{iwire}(trl,:) = bSpikeeKF{1,iwire}(trl,:)
        end
        else
        fprintf('No wave-clus for this wire.\n');
        end
    end
    
end

%save ('Rastertime_eachtrial_posc1', 'eR', 'eKF')
%save ('Rastertime_eachtrial_negc1', 'eR', 'eKF')
%save ('Rastertime_eachtrial_posc2', 'eR', 'eKF')
%save ('Rastertime_eachtrial_negc2', 'eR', 'eKF')
%%
clear all
close all
clc

bins = '20_bins'

cd (['C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\Prestim\20_bins']);
load ('Rastertime_eachtrial_posc1')
posc1eR = eR
posc1eKF = eKF

clear eR
clear eKF

load ('Rastertime_eachtrial_negc1')
negc1eR = eR
negc1eKF = eKF

clear eR
clear eKF

load ('Rastertime_eachtrial_posc2')
posc2eR = eR
posc2eKF = eKF

clear eR
clear eKF

load ('Rastertime_eachtrial_negc2')
negc2eR = eR
negc2eKF = eKF

for s=1:7
    AlleRrastertime{:,s} = [posc1eR{1,s} negc1eR{1,s} posc2eR{1,s}]
    AlleKFrastertime{:,s} = [posc1eKF{1,s} negc1eKF{1,s} posc2eKF{1,s}]
end       

A= AlleRrastertime
B= AlleKFrastertime
        
for s = 1:7
    for iw= 1:size(A{s},2)
    for trl = 1:size(A{s}{1,iw},1)
    Asum{s}{1,iw}(trl,:) = sum(A{s}{1,iw}(trl,:))
    sumiwi{s}(1,iw) = sum(Asum{s}{1,iw})
    end
    end
end

for s=1:7
[I,J,V] = find(sumiwi{s} ~= 0)
tokeep{s} = J
end

As1 = A{1}(1,tokeep{1});
As2 = A{2}(1,tokeep{2});
As3 = A{3}(1,tokeep{3});
As4 = A{4}(1,tokeep{4});
As5 = A{5}(1,tokeep{5});
As6 = A{6}(1,tokeep{6});
As7 = A{7}(1,tokeep{7});

Bs1 = B{1}(1,tokeep{1});
Bs2 = B{2}(1,tokeep{2});
Bs3 = B{3}(1,tokeep{3});
Bs4 = B{4}(1,tokeep{4});
Bs5 = B{5}(1,tokeep{5});
Bs6 = B{6}(1,tokeep{6});
Bs7 = B{7}(1,tokeep{7});



Au = [As1 As2 As3 As4 As5 As6 As7]
Bu = [Bs1 Bs2 Bs3 Bs4 Bs5 Bs6 Bs7]

%save('allun32_rastertimeprestim_trl.mat', 'Au','Bu'); 
%%
clear all
close all
clc

cd C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\Prestim\20_bins
load allun32_rastertimeprestim_trl.mat
load time

toi = [-0.5 1.5];

pt1 = nearest(LFP_res.time{1},toi(1));
pt2 = nearest(LFP_res.time{1},toi(2));

for u=1:32
Acount = Au{u}(:,[pt1:pt2]);
totA(u) = sum(sum(Acount,2),1);
end

%%
[c_sorted, c_order] = sort(totA,'descend');

Au = Au(:,c_order);
Bu = Bu(:,c_order);

%%
for u=1:32
clear Acount
Acount = Au{u};
totAtime(u,:) = mean(Acount,1);
end

for u=1:32
clear Bcount
Bcount = Bu{u};
totBtime(u,:) = mean(Bcount,1);
end
%%
NbspikeeR = totAtime([1:31],[pt1:pt2]);
NbspikeeKF = totBtime([1:31],[pt1:pt2]);

%% sum of spikes /count
tb= [0:160:4000]
tbins= tb(2:end)

for u=1:31
    NbspikeA = NbspikeeR(u,:)
    NbspikeA = NbspikeA'
    v=0
    for i=tbins
        v=v+1
        clear var
        if v==1
            var = NbspikeA(v:tbins(v))
            NbspikeA_bin = sum(var,1)
            A(v,:)= NbspikeA_bin
        else
            var = NbspikeA(tbins(v-1):tbins(v))
            NbspikeA_bin = sum(var,1)
            A(v,:)= NbspikeA_bin
        end
    end
   CountA(u,:) = A
end

%% sum of spikes /count
tb= [0:160:4000]

tbins= tb(2:end)

for u=1:31
    NbspikeB = NbspikeeKF(u,:)
    NbspikeB = NbspikeB'
    v=0
    for i=tbins
        v=v+1
        clear var
        if v==1
            var = NbspikeB(v:tbins(v))
            NbspikeB_bin = sum(var,1)
            B(v,:)= NbspikeB_bin
        else
            var = NbspikeB(tbins(v-1):tbins(v))
            NbspikeB_bin = sum(var,1)
            B(v,:)= NbspikeB_bin
        end
    end
   CountB(u,:) = B
end
%% supplementary Fig.22
figure
meanA = mean(CountA)
meanB = mean(CountB)
FRA = meanA*1000/80;% 25 bins
FRB = meanB*1000/80
subplot(1,2,1)
bar(FRA,'FaceColor', [1, 0, 0], 'EdgeColor', 'none');
set(gca,'XTick',1:1.25:25,'XTickLabel',-0.5:0.1:1.5);
xlabel('Time (s)');
ylabel ('mean FR (Hz)')
hold on
subplot(1,2,2)
bar(FRB,'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
set(gca,'XTick',1:1.25:25,'XTickLabel',-0.5:0.1:1.5);
xlabel('Time (s)');
ylabel ('mean FR (Hz)')
hold on
%% raster plot
%% plot er and ekf
fig = figure('rend','painters','pos',[10 10 400 700])
%%
toir = [-0.5 1.5];
pt1r = nearest(LFP_res.time{1},toir(1));
pt2r = nearest(LFP_res.time{1},toir(2));
timer= LFP_res.time{1}(1,[pt1r:pt2r])
%%
toig = [0.4102 1.1006];
pt1g = nearest(timer,toig(1));% in a vector from 1:8000
pt2g = nearest(timer,toig(2));
%%
for u=1:8;
    Autoplot = Au{u};
    Butoplot = Bu{u};
    
    binarySpikes = []
    binarySpikes = Autoplot;
    
    h_raster = subplot(8,1,u)
    position = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
    pos = position(u)
    set(h_raster, 'Position', [0.05, pos, 0.32, 0.08]);
    MarkerFormat =  struct()
    MarkerFormat.Color = [1 0 0]
    plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
    hold on
end

%%
fig = figure('rend','painters','pos',[10 10 400 700])

for u=1:8;
Autoplot = Au{u};
Butoplot = Bu{u};

binarySpikes = []
binarySpikes = Butoplot;
h_raster1 = subplot(8,1,u)
position1 = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
pos1 = position1(u)
set(h_raster1, 'Position', [0.52, pos1, 0.32, 0.08]);
MarkerFormat =  struct()
MarkerFormat.Color = [1 0 1]
plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
hold on
end

%%
fig = figure('rend','painters','pos',[10 10 400 700])

for u=9:16;
Autoplot = Au{u};
Butoplot = Bu{u};

binarySpikes = []
binarySpikes = Autoplot;
t=u-8
h_raster1 = subplot(8,1,t)
position1 = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
pos1 = position1(t)
set(h_raster1, 'Position', [0.05, pos1, 0.32, 0.08]);
MarkerFormat =  struct()
MarkerFormat.Color = [1 0 0]
plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
hold on
end
%%

fig = figure('rend','painters','pos',[10 10 400 700])

for u=9:16;
Autoplot = Au{u};
Butoplot = Bu{u};

binarySpikes = []
binarySpikes = Butoplot;
t=u-8
h_raster1 = subplot(8,1,t)
position1 = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
pos1 = position1(t)
set(h_raster1, 'Position', [0.52, pos1, 0.32, 0.08]);
MarkerFormat =  struct()
MarkerFormat.Color = [1 0 1]
plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
hold on
end

%%
fig = figure('rend','painters','pos',[10 10 400 700])

for u=17:24;
Autoplot = Au{u};
Butoplot = Bu{u};

binarySpikes = []
binarySpikes = Autoplot;
t=u-16
h_raster1 = subplot(8,1,t)
position1 = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
pos1 = position1(t)
set(h_raster1, 'Position', [0.05, pos1, 0.32, 0.08]);
MarkerFormat =  struct()
MarkerFormat.Color = [1 0 0]
plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
hold on
end
%%

fig = figure('rend','painters','pos',[10 10 400 700])

for u=17:24;
Autoplot = Au{u};
Butoplot = Bu{u};

binarySpikes = []
binarySpikes = Butoplot;
t=u-16
h_raster1 = subplot(8,1,t)
position1 = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
pos1 = position1(t)
set(h_raster1, 'Position', [0.52, pos1, 0.32, 0.08]);
MarkerFormat =  struct()
MarkerFormat.Color = [1 0 1]
plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
hold on
end

%%
fig = figure('rend','painters','pos',[10 10 400 700])

for u=25:32;
Autoplot = Au{u};
Butoplot = Bu{u};

binarySpikes = []
binarySpikes = Autoplot;
t=u-24
h_raster1 = subplot(8,1,t)
position1 = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
pos1 = position1(t)
set(h_raster1, 'Position', [0.05, pos1, 0.32, 0.08]);
MarkerFormat =  struct()
MarkerFormat.Color = [1 0 0]
plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
hold on
end
%%

fig = figure('rend','painters','pos',[10 10 400 700])

for u=25:32
Autoplot = Au{u};
Butoplot = Bu{u};

binarySpikes = []
binarySpikes = Butoplot;
t=u-24
h_raster1 = subplot(8,1,t)
position1 = [0.84 0.74 0.64 0.54 0.44 0.34 0.24 0.14]
pos1 = position1(t)
set(h_raster1, 'Position', [0.52, pos1, 0.32, 0.08]);
MarkerFormat =  struct()
MarkerFormat.Color = [1 0 1]
plotSpikeRaster(binarySpikes(:,[pt1:pt2]),'PlotType','scatter','XLimForCell',[0 0.201],'MarkerFormat',MarkerFormat);
    %title(['eR hc spikes', ' unit ',num2str(u)]);
    set(gca,'XTick',1:1000:4000,'XTickLabel',-0.5:0.5:1.5);
    %xlabel('Time (ms)');
    %ylabel('Trial');
    set(gca,'xtick',[])
    hold on
    y = ylim; % current y-axis limits
    plot([pt1g pt1g],[y(1) y(2)],'k-')
    hold on
    plot([pt2g pt2g],[y(1) y(2)],'k-')
hold on
end

