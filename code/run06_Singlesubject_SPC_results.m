
clear all; close all; clc;

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','circstat-matlab'))
addpath(fullfile(oripath,'Costalozanoetal','code','pacoi'))

direction='posSpikes';% posSpikes | negSpikes 
cluster= 'Cluster1';% Cluster1 | Cluster2


bins = '20_bins' % 

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

angle_bins      = deg2rad(-180:18:180);
angle_centers  	= transpose(movmean(angle_bins, 2, 'endpoints', 'discard'));


v=0

iclus = 1% change this value when you want to look at cluster 2
for isub = 1:length(subjects)
    
    freq_band = [3,12];
    %for 
        fb = 1
        v=v+1;
        cd (['C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\' cluster '\' direction]);
        name = ['freq',num2str(freq_band(fb,1)),'to',num2str(freq_band(fb,2))]
        cd (name)
        load(['Z',subjects{isub},'_data']);
        
        % AlleR contain all spikes occurring in this condition for each
        % wire, some are empty because no spikes were found
        for iwire = 1:size(SPCeR,2)
            cnt=0
            for trl= 1:size(SPCeR{iwire},2)
                for i= 1:size(SPCeR{iwire}{trl},1)
                    cnt= cnt+1
                    AlleR{isub}{iwire}(cnt) = SPCeR{iwire}{trl}(i,:)
                end
            end
         if cnt == 0
         AlleR{isub}{iwire} = []
         end
        end
        
        clear count 
        clear iwire
         for iwire = 1:size(SPCeKF,2)
            cnt=0
            for trl= 1:size(SPCeKF{iwire},2)
                for i= 1:size(SPCeKF{iwire}{trl},1)
                    cnt= cnt+1
                    AlleKF{isub}{iwire}(cnt) = SPCeKF{iwire}{trl}(i,:)
                end
            end
         if cnt == 0
         AlleKF{isub}{iwire} = []
         end
         end
        
          toi = [0.4102 1.1006];
          t1 = round(toi(1).*2000/1); % before I run from 0 to 1 
          t2 = round(toi(2).*2000/1);
         
         % calculate how many spikes occured in each condition in a window
         % from 0.41 to 1.1 sec 
         for iw = 1:size(bSpikeeR,2);
                for trl = 1:size(bSpikeeR{iw},1)
                    [i,j,s] = find(bSpikeeR{iw}(trl,[t1:t2]) == 1);
                    NbSpikeseR{isub}(trl,iw) = sum(s);
                end
            end
            
            for iw = 1:size(bSpikeeKF,2);
                for trl = 1:size(bSpikeeKF{iw},1)
                    [i2,j2,s2] = find(bSpikeeKF{iw}(trl,[t1:t2]) == 1);
                    NbSpikeseKF{isub}(trl,iw) = sum(s2);
                end
            end
         
         
         for iwire = 1 :size(AlleR{isub},2)
             
             fig1= figure;
             subplot(1,2,1)
             polarhistogram(AlleR{isub}{iwire},20,'Normalization','probability','FaceColor','red','FaceAlpha',.5);
             title ('eR');
             hold on
             subplot(1,2,2)
             polarhistogram(AlleKF{isub}{iwire},20,'Normalization','probability','FaceColor','black','FaceAlpha',.5);
             title ('eKF');
             axes( 'Position', [0, 0.1, 1, 0.045]) ;
             set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
             text(0.5, 0,['n eR = ', num2str(size(AlleR{1,isub}{1,iwire},2)), ' / n eKF = ', num2str(size(AlleKF{1,isub}{1,iwire},2))], ...
                 'fontunits', 'centimeters', 'fontsize', 0.3, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom'); % number of spikes
             axes( 'Position', [0, 0.95, 1, 0.05] ) ;
             set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
             text( 0.5, 0, [subjects{isub},' - SPC Amy theta (3-12 Hz) Hc - wire' num2str(iwire),' - ', direction,' - ', ' cluster' num2str(iclus)], 'FontSize', 12, 'FontWeight', 'Bold', ...
                 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
             
            
             
             % btpheR contain nb of spikes occuring in each of the 20 bins
             % for all trial of one condition
             if size(binnedThetaPhaseseR{iwire},2) ~= 0
                 for rp = 1:size(binnedThetaPhaseseR{iwire},2);
                     btpheR(:,rp)= binnedThetaPhaseseR{1,iwire}{1,rp}
                 end
             else
                 btpheR = []
             end
            
            if size(binnedThetaPhaseseKF{iwire},2) ~= 0
            for rp = 1:size(binnedThetaPhaseseKF{iwire},2);
            btpheKF(:,rp)= binnedThetaPhaseseKF{1,iwire}{1,rp}  
            end
            else
                 btpheKF = []
            end
            
            % sum over trials spikes for each bin
            btpheRall = sum(btpheR,2);
            btpheKFall = sum(btpheKF,2);

            
            %% polar histogram to depict theta modulation
            if size(btpheRall,1) ~= 0
            fig2= figure 
            set(gcf, 'Position', get(0, 'Screensize'));
            subplot(1,3,1)
            x = []; y= []
            x_eR = []; y_eR = []
            radii = 0.25:0.25:1;
            for pxd = 1:numel(radii)
                plot(cos(0:0.001:2*pi) .* radii(pxd), sin(0:0.001:2*pi) .* radii(pxd), 'k--')
            end
            [x, y] = pol2cart(angle_centers, btpheRall ./ max(btpheRall));
            patch([x; x(1)], [y; y(1)], [1 0 0], ...
                'facealpha', 0.2, ...
                'edgecolor', [1 0 0], ...
                'linewidth', 1);
            % circ mean of eR over bins
            [x_eR, y_eR] = pol2cart(circ_mean(angle_centers, btpheRall), 1);
            hold on
            plot([0, x_eR], [0, y_eR], 'r-', ...
                'linewidth', 3.5);
            text(0.75, -0.85, ['n eR = ', num2str(size(AlleR{1,isub}{1,iwire},2))], ...
                'fontunits', 'centimeters', 'fontsize', 0.3); % number of spikes
            text(1.1, 0, '0°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'left', 'verticalalignment', 'middle'); % 0° label
            text(-1.1, 0, [0177, '180°'], ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'right', 'verticalalignment', 'middle'); % 180° label
            text(0, 1.1, '90°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'center', 'verticalalignment', 'bottom'); % 90° label
            text(0, -1.1, '-90°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'center', 'verticalalignment', 'top'); % 90° label
            set(gca, ...
                'xlim', [-1, 1], 'ylim', [-1, 1]);
            axis square; axis off;
            hold on
            subplot(1,3,2)
            x1 = []; y1= []
            x_eKF = []; y_eKF = []
            radii = 0.25:0.25:1;
            for pxd = 1:numel(radii)
                plot(cos(0:0.001:2*pi) .* radii(pxd), sin(0:0.001:2*pi) .* radii(pxd), 'k--')
            end
            [x1, y1] = pol2cart(angle_centers, btpheKFall ./ max(btpheKFall));
            patch([x1; x1(1)], [y1; y1(1)], [0.1 0.1 0.1], ...
                'facealpha', 0.2, ...
                'edgecolor', [0.1 0.1 0.1], ...
                'linewidth', 1);
            % circ mean eKF over bins
            [x_eKF, y_eKF] = pol2cart(circ_mean(angle_centers, btpheKFall), 1);
            hold on
            plot([0, x_eKF], [0, y_eKF], 'k-', ...
                'linewidth', 3.5);
            text(0.75, -0.85, ['n eKF = ', num2str(size(AlleKF{1,isub}{1,iwire},2))], ...
                'fontunits', 'centimeters', 'fontsize', 0.3); % number of spikes
            text(1.1, 0, '0°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'left', 'verticalalignment', 'middle'); % 0° label
            text(-1.1, 0, [0177, '180°'], ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'right', 'verticalalignment', 'middle'); % 180° label
            text(0, 1.1, '90°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'center', 'verticalalignment', 'bottom'); % 90° label
            text(0, -1.1, '-90°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'center', 'verticalalignment', 'top'); % 90° label
            set(gca, ...
                'xlim', [-1, 1], 'ylim', [-1, 1]);
            axis square; axis off;
            hold on
            subplot(1,3,3)
            radii = 0.25:0.25:1;
            for pxd = 1:numel(radii)
                plot(cos(0:0.001:2*pi) .* radii(pxd), sin(0:0.001:2*pi) .* radii(pxd), 'k--')
            end
            hold on
            % realign the two conditions
            ang_bineR = []; ang_bineKF = []
            x_eRcent = []; y_eRcent = []
            x_eKFcent = []; y_eKFcent = []
            ang_bineR = btpheRall'
            ang_bineKF =  btpheKFall'
            x_dbin =angle_centers'
            [angbin1a, angbin2a2, flag(s),angbin2a] = realign2phases(ang_bineR,ang_bineKF,x_dbin);
            eMa{isub}(iwire,:) = angbin1a;
            nMa{isub}(iwire,:) = angbin2a2;
  
            [x_eRcent, y_eRcent] = pol2cart(circ_mean(angle_centers, angbin1a'), 1);
            [x_eKFcent, y_eKFcent] = pol2cart(circ_mean(angle_centers, angbin2a2'), 1);
            hold on
            plot([0, x_eRcent], [0, y_eRcent], 'r-', ...
                'linewidth', 3.5);
            hold on
            plot([0, x_eKFcent], [0, y_eKFcent], 'k-', ...
                'linewidth', 3.5);
            text(1.1, 0, '0°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'left', 'verticalalignment', 'middle'); % 0° label
            text(-1.1, 0, [0177, '180°'], ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'right', 'verticalalignment', 'middle'); % 180° label
            text(0, 1.1, '90°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'center', 'verticalalignment', 'bottom'); % 90° label
            text(0, -1.1, '-90°', ...
                'fontunits', 'centimeters', 'fontsize', 0.3, ...
                'horizontalalignment', 'center', 'verticalalignment', 'top'); % 90° label
            set(gca, ...
                'xlim', [-1, 1], 'ylim', [-1, 1]);
            axis square; axis off;
            hold on
             axes( 'Position', [0, 0.95, 1, 0.05] ) ;
             set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%              text(0, -1.5,[num2str(real(zN_all(u)),'%3.2f'), ' - ', num2str(real(zE_all(u)),'%3.2f'), ' = ', num2str(real(popptt(u)),'%3.2f'), ' rad'], ...
%         'fontunits', 'centimeters', 'fontsize', 0.3, 'fontweight', 'bold' ,'HorizontalAlignment', 'Center','VerticalAlignment','bottom');
         
             text( 0.5, 0, [subjects{isub},' - SPC Amy theta (3-12 Hz) Hc - wire' num2str(iwire),' - ', direction,' - ', ' cluster' num2str(iclus)], 'FontSize', 12, 'FontWeight', 'Bold', ...
                 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
             
            %saveas(fig2,['Z',subjects{isub},'_SPCconds_meancent',num2str(iwire)],'png')
            %saveas(fig2,['Z',subjects{isub},'_SPCconds_meancent',num2str(iwire)],'fig')
            
            else
            fprintf('No wave-clus for this wire.\n');
            end
            
                        
  clear btpheR; clear btpheKF; clear btpheRall; clear btpheKFall; clear angbin1a; clear angbin2a2;
  
  %file = ['Phase_realign',direction,cluster,'.mat'];
  %save(file,'eMa','nMa');   % these are eR and eKF after phase realignment      
  
 close all
 cd 
       end
end
%% prepare group results
clear all
close all
clc

bins = '20_bins' ;

% to account for all units I load eR and eKF (phase realigned) for each subject and
% for cluster 1 pos spike

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','SingleUnit','Phase_realignposSpikesCluster1.mat'))

for s=1:size(eMa,2)
    iwire = 1:size(eMa{s},1)
    sumiwi{s}(iwire) = sum(eMa{s}(iwire,:),2); 
    [I,J,V] = find(sumiwi{s} ~= 0)
    tokeep{s} = J
end

% collapse results over subjects and keep only wires with action potential 
AlleMa = [eMa{1}(tokeep{1},:); eMa{2}(tokeep{2},:); eMa{3}(tokeep{3},:); eMa{4}(tokeep{4},:);eMa{5}(tokeep{5},:); eMa{6}(tokeep{6},:); eMa{7}(tokeep{7},:)]
AllnMa = [nMa{1}(tokeep{1},:); nMa{2}(tokeep{2},:); nMa{3}(tokeep{3},:); nMa{4}(tokeep{4},:);nMa{5}(tokeep{5},:); nMa{6}(tokeep{6},:); nMa{7}(tokeep{7},:)]

eMa1 = AlleMa;
nMa1 = AllnMa;

%save ('pos_c1.mat','eMa1','nMa1');
%save ('utokeep_pc1.mat','tokeep')
%%
clear all
close all
clc

% for cluster 1 neg spike
oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','SingleUnit','Phase_realignnegSpikesCluster1.mat'))


for s=1:size(eMa,2)
    iwire = 1:size(eMa{s},1)
    sumiwi{s}(iwire) = sum(eMa{s}(iwire,:),2); 
    [I,J,V] = find(sumiwi{s} ~= 0)
    tokeep{s} = J
end

AlleMa = [eMa{1}(tokeep{1},:); eMa{2}(tokeep{2},:); eMa{3}(tokeep{3},:); eMa{4}(tokeep{4},:);eMa{5}(tokeep{5},:); eMa{6}(tokeep{6},:); eMa{7}(tokeep{7},:)]
AllnMa = [nMa{1}(tokeep{1},:); nMa{2}(tokeep{2},:); nMa{3}(tokeep{3},:); nMa{4}(tokeep{4},:);nMa{5}(tokeep{5},:); nMa{6}(tokeep{6},:); nMa{7}(tokeep{7},:)]


eMa2 = AlleMa;
nMa2 = AllnMa;

%save ('neg_c1.mat','eMa2','nMa2');
%save ('utokeep_nc1.mat','tokeep')

%%
clear all
close all
clc

% for cluster 2 pos spike
oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','SingleUnit','Phase_realignposSpikesCluster2.mat'))

% find wire that do not have any action potential
for s=1:size(eMa,2)
    iwire = 1:size(eMa{s},1)
    sumiwi{s}(iwire) = sum(eMa{s}(iwire,:),2); 
    [I,J,V] = find(sumiwi{s} ~= 0)
    tokeep{s} = J
end

% erase empty wires and put together results for all subjects
AlleMa = [eMa{1}(tokeep{1},:); eMa{2}(tokeep{2},:); eMa{3}(tokeep{3},:); eMa{4}(tokeep{4},:);eMa{5}(tokeep{5},:); eMa{6}(tokeep{6},:); eMa{7}(tokeep{7},:)]
AllnMa = [nMa{1}(tokeep{1},:); nMa{2}(tokeep{2},:); nMa{3}(tokeep{3},:); nMa{4}(tokeep{4},:);nMa{5}(tokeep{5},:); nMa{6}(tokeep{6},:); nMa{7}(tokeep{7},:)]


eMa3 = AlleMa;
nMa3 = AllnMa;


%save ('pos_c2.mat','eMa3','nMa3');
%save ('utokeep_pc2.mat','tokeep')

%%
clear all
close all
clc

% for cluster 2 neg spike
oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
load(fullfile(oripath,'data2github','SingleUnit','Phase_realignnegSpikesCluster2.mat'))

for s=1:size(eMa,2)
    iwire = 1:size(eMa{s},1)
    sumiwi{s}(iwire) = sum(eMa{s}(iwire,:),2); 
    [I,J,V] = find(sumiwi{s} ~= 0)
    tokeep{s} = J
end

AlleMa = [eMa{1}(tokeep{1},:); eMa{2}(tokeep{2},:); eMa{3}(tokeep{3},:); eMa{4}(tokeep{4},:);eMa{5}(tokeep{5},:); eMa{6}(tokeep{6},:); eMa{7}(tokeep{7},:)]
AllnMa = [nMa{1}(tokeep{1},:); nMa{2}(tokeep{2},:); nMa{3}(tokeep{3},:); nMa{4}(tokeep{4},:);nMa{5}(tokeep{5},:); nMa{6}(tokeep{6},:); nMa{7}(tokeep{7},:)]


eMa4 = AlleMa;
nMa4 = AllnMa;

% save ('neg_c2.mat','eMa4','nMa4'); is empty