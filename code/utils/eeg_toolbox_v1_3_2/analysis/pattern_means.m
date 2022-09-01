function eeg = pattern_means(eeg, resDir, regressors, interactions, overall)

if ~exist(fullfile(resDir, 'data'), 'dir')
  mkdir(fullfile(resDir, 'data'));
end

eeg.resDir = resDir;

for s=1:length(eeg.subj)
  fprintf('\n%s\n', eeg.subj(s).id);
  
  eeg.subj(s).patmeansFile = fullfile(resDir, 'data', [eeg.subj(s).id '_patmeans.mat']);
  
  if ~lockFile(eeg.subj(s).patmeansFile)
    save(fullfile(eeg.resDir, 'eeg.mat'), 'eeg');
    continue
  end
  
  load(eeg.subj(s).patFile);
  load(eeg.subj(s).regFile);
  
  reg = filterStruct(reg, 'ismember(name, varargin{1})', regressors);
  
  % get mean values for each regressor
  for i=1:length(reg)
    patmeans(i).name = reg(i).name;
    patmeans(i).vals = unique(reg(i).vec);    
    for j=1:length(patmeans(i).vals)
      patmeans(i).mat{j} = NaN(length(eeg.subj(s).chan), size(pat,1), size(pat,2));
    end
  end
  
  if interactions
    intermeans(1).name = [patmeans(1).name ' X ' patmeans(2).name];
    for i=1:length(patmeans(1).vals)
      for j=1:length(patmeans(2).vals)
	intermeans(1).vals{i,j} = [patmeans(1).vals(i) patmeans(2).vals(j)];
	intermeans(1).mat{i,j} = NaN(length(eeg.subj(s).chan), size(pat,1), size(pat,2));
      end
    end
  end
  
  for c=1:length(eeg.subj(s).chan)
    fprintf('%d ', eeg.subj(s).chan(c).number);
  
    for f=1:size(pat, 1)
      for b=1:size(pat, 2)
	thispat = pat(f,b).mat(:,c);
	
	for i=1:length(reg)
	  for j=1:length(patmeans(i).vals)
	    if iscell(patmeans(i).vals)
	      thiscond = thispat(strcmp(reg(i).vec, patmeans(i).vals{j}));
	    else
	      thiscond = thispat(reg(i).vec==patmeans(i).vals(j));
	    end
	    patmeans(i).mat{j}(c,f,b) = nanmean(thiscond);
	  end
	end
	
	if overall
	  patmeans(length(reg)+1).name = 'overall';
	  patmeans(length(reg)+1).vals = [];
	  patmeans(length(reg)+1).mat{1}(c,f,b) = nanmean(thispat);
	end
	
	if interactions
	  for i=1:length(patmeans(1).vals)
	    for j=1:length(patmeans(2).vals)
	      thiscond = reg(1).vec==patmeans(1).vals(i) & reg(2).vec==patmeans(2).vals(j);
	      intermeans(1).mat{i,j}(c,f,b) = nanmean(thispat(thiscond));
	    end
	  end
	end
	
      end
    end
  end
  
  for i=1:length(reg)
    for j=1:length(patmeans(i).vals)
      patmeans(i).mat{j} = squeeze(patmeans(i).mat{j});
    end
  end
  
  save(eeg.subj(s).patmeansFile, 'patmeans', 'intermeans');
  releaseFile(eeg.subj(s).patmeansFile);
  save(fullfile(eeg.resDir, 'eeg.mat'), 'eeg');
  
end