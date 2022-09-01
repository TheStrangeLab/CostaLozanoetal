function [clusPvalP,clusP,clusPvalN,clusN] = getSignifClusters(p_sig,p_boot,pthresh)
% FUNCTION getSignifClusters - find the significant clusters within
% the data (based on Maris & Oostenveld (2007) J Neurosci Meth
% 164:177-190)
%
% function [clusPvalP,clusP,clusPvalN,clusN] =
% getSignifClusters(p_sig,p_boot,pthresh)
%
% INPUT ARGS:
% psig = significance matrix (channels,frequency,time)
% pboot = matrix of bootstrapped p-values
% (Nperm,chan,frequency,time) where Nperm is the number of
% bootstrap iterations
% pthresh=0.05  - significance threshold at which to threshold the data
%
% OUTPUT ARGS
% structure with significant clusters
%
% this function uses the findCluster.m function from Fieldtrip,
% which relies on the Image Processing toolbox
%

if ~exist('pthresh','var')
  thresh = 0.05;
end


% threshold the data
onoffSpos = p_sig<pthresh;
onoffSneg = p_sig>(1-pthresh);

% prepare the neighborhood structure
%% for sEEG
cfg = [];
[cfg.elec,label] = prepChanStruct;
cfg.neighbourdist = 3;
neighbors = neighbourselection(cfg);
neighborMat = makechanneighbstructmat(neighbors,label);

[clusPos,numPos]=findcluster(onoffSpos,neighborMat);
[clusNeg,numNeg]=findcluster(onoffSneg,neighborMat);
% get the summed cluster stats for every cluster
clusStatsP = zeros(1,numPos); % clus = chan x freq x time
clusStatsN = zeros(1,numNeg);
% the cluster statistic is the sum of the norminv p-values
for c = 1:numPos
  clusStatsP(c) = sum(norminv(p_sig(find(clusPos==c))));
end
for c = 1:numNeg
  clusStatsN(c) = sum(norminv(p_sig(find(clusNeg==c))));
end

% get the clusters for every bootstrap iteration
iterations = size(p_boot,1);
clusCollectP = zeros(1,iterations);
clusCollectN = zeros(1,iterations);
for i = 1:size(p_boot,1)
  fprintf('%d ',i);
  onoffBpos = squeeze(p_boot(i,:,:,:))<pthresh;
  [clusBpos,numBpos] = findcluster(onoffBpos,neighborMat);
  onoffBneg = squeeze(p_boot(i,:,:,:))>(1-pthresh);
  [clusBneg,numBneg] = findcluster(onoffBneg,neighborMat);
  bootClusStatP = zeros(1,numBpos);
  bootClusStatN = zeros(1,numBneg);
  % get the summed cluster stats
  for c = 1:numBpos
    bootClusStatP(c) = sum(norminv(p_boot(i,find(clusBpos==c))));
  end
  for c = 1:numBneg
    bootClusStatN(c) = sum(norminv(p_boot(i,find(clusBneg==c))));
  end
  % get the maximum cluster stat for this iteration
  clusCollectP(i) = min(bootClusStatP);
  clusCollectN(i) = max(bootClusStatN);
end

clusCollectP = sort(clusCollectP);
clusCollectN = sort(clusCollectN);
% now compare the obtained clusters to the bootstrapped clusters
clusPvalP = ones(numPos,1);
clusPvalN = ones(numNeg,1);
for c = 1:numPos
  % check the p-values in both positive and negative directions
  pvalP = 1; pvalN = 1;
  d = find(clusStatsP(c)<clusCollectP);
  if ~isempty(d)
    clusPvalP(c) = (d(1)/iterations);
  end
end

for c = 1:numNeg
  d = find(clusStatsN(c)>clusCollectN);
  if ~isempty(d)
    clusPvalN(c) = 1-(d(end)/iterations);
  end  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN SUBFUNCTION MAKECHANNEIGHBSTRUCTMAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [channeighbstructmat] = makechanneighbstructmat(neighbours,label);

% MAKECHANNEIGHBSTRUCTMAT makes the makes the matrix containing the channel
% neighbourhood structure.

nchan=length(label);
channeighbstructmat = logical(zeros(nchan,nchan));
if length(neighbours)==0
    error('Empty neighborhood structure ''cfg.neighbours''.');
end;
for chan=1:length(neighbours)
    [seld] = match_str(label, neighbours{chan}.label);
    nblabel = neighbours{chan}.neighblabel;
    nblabel = [nblabel{:}];
    [seln] = match_str(label, nblabel);
    if isempty(seld)
        % this channel was not present in the data
        continue;
    else
        % add the neighbours of this channel to the matrix
        channeighbstructmat(seld, seln) = 1;
    end
end;
