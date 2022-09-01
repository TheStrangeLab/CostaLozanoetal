function data = prepare_iapsKF_9patients(data)
channels = cell(1,9);
channels{1} = {'A1';'A2';'HC1'};
channels{3} = {'A1';'A2';'HC1'};
channels{9} = {'A1';'A2';'HC1'};
channels{2} = {'A1';'A2';'HC1';'HC2'};
channels{4} = {'A1';'HC1';'HC2'};
channels{5} = {'A1';'HC1';'HC2'};
channels{6} = {'A1';'HC1'};
channels{7} = {'A1';'HC1'};
channels{8} = {'A1';'HC1'};

for s=1:9
  data{s}.label = channels{s};
  if size(data{s}.label,1)~=size(data{s}.trial{1},1)
    error('wrong channel number')
  end
end