function [trl flag nrpt triggerpoints data fsample] = mytrialfun_IntraCranial_iaps(cfg)
%
hdr = ft_read_header(cfg.dataset);
fsample = hdr.Fs;
data = ft_read_data(cfg.dataset,'chanindx',1); %read only "Event" channel
% data = reshape(data,size(data,1),size(data,2)*size(data,3)); %reshape data to continuous data representation

% conn = bwconncomp(data);
% inds = conn.PixelIdxList;
%
% triggerpoints = zeros(size(inds));
% for k=1:size(inds,2)
%     triggerpoints(1,k) = inds{k}(1,1);% select the first time sample
% end

triggerpoints = find(diff([data(1) data(:)'])>0); %get trigger flanks (up)

% 21-Aug-2020 01:35:40 - DLS
% check whether trigger points have more than 1 sample. In some patients
% the tigger was on 2-3 samples
iti = diff(triggerpoints);
consmp = find(iti < 5); % 5 samples as tolerance
torest = iti(consmp)-1; % values to subtract

todelete = [];
for k=1:size(torest,2)
  todelete = [todelete consmp(k)+torest(k)];
end

triggerpoints(todelete)=[];

if size(triggerpoints,2)==size(cfg.triggers,2)
  trl = [];
  flag = 1;
  for i=1:length(triggerpoints)
    begsample = triggerpoints(i) - cfg.prestim*hdr.Fs;
    endsample = triggerpoints(i) + cfg.poststim*hdr.Fs - 1;
    offset = -cfg.prestim*hdr.Fs;
    trl(end+1, :) = round([begsample endsample offset]);
  end
  trl(:,end+1) = cfg.triggers';
  nrpt = size(trl,1);
else
%     error('mismatch between EDF and number of triggers')
  flag = 0;
  trl=[];
  nrpt = size(triggerpoints,2);
end
