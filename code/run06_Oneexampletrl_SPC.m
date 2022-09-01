% This script reproduces an example trial showing hippocampal spikes and amygdala theta phase as in Supplementary Fig. 23
clear all; close all; clc;

% add fieldtrip
restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','eeg_toolbox_v1_3_2'))% kahana eeg toolbox
addpath(fullfile(oripath,'Costalozanoetal','code','utils','MatlabImportExport_v6.0.0'))% 
warning('off', 'all'); % suppress all warning

%% t1 isub=4 trl 13
subjects    = {...
    'Patient1'; ...
    'Patient2'; ...
    'Patient200'; ...
    'Patient5'; ...
    'Patient6'; ...
    'Patient8'; ...
    'Patient10'; ...
    };

isub=4

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','SingleUnit','Prestim',['Z',subjects{isub},'_data_exampletrl'])) 

toi = [0.001 3.2];
t1 = round(toi(1).*2000/1); 
t2 = round(toi(2).*2000/1);

figure
iwire = 1
trl = 13;
plot(bSpikeeR{iwire}(trl,[t1:t2]),'b','LineWidth', 1)
hold on
plot(bandphaseseR{trl}(1,[t1:t2]),'r')
ylabel ('amygdala theta phase')
xlabel ('time in sec')
set(gca,'XTick',1:400:6400,'XTickLabel',0:0.2:3.2);
set(gcf,'color','w');
box off
% add legend
Lgnd = legend ('hippocampal spikes','theta phase (3-12 hz)')
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.01;

hold on
plot(bSpikeeKF{iwire}(trl,[t1:t2]),'b','LineWidth', 1)
hold on
plot(bandphaseseKF{trl}(1,[t1:t2]),'m')
ylabel ('amygdala theta phase')
xlabel ('time in sec')
set(gca,'XTick',1:400:6400,'XTickLabel',0:0.2:3.2);
set(gcf,'color','w');
box off
% add legend
Lgnd = legend ('hippocampal spikes','theta phase (3-12 hz)')
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.01;
%%
figure
for i=1:15
subplot(5,3,i)
iwire = 1
trl = 13;
plot(bSpikeeR{iwire}(i,[t1:t2]),'b','LineWidth', 1)
hold on
plot(bandphaseseR{i}(1,[t1:t2]),'r')
ylabel ('amygdala theta phase')
xlabel ('time in sec')
set(gca,'XTick',1:400:6400,'XTickLabel',0:0.2:3.2);
set(gcf,'color','w');
box off
% add legend
Lgnd = legend ('hippocampal spikes','theta phase (3-12 hz)')
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.01;
end