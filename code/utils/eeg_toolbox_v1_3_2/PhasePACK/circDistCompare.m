function p=circDistCompare2(varargin)
%Nonparametric test to compare circular distributions and indicate if
%they're identical.
%
%Taken from Fisher 5.3.6
%coded by JJ
%Note: this hasn't been extensively tested yet....

if length(varargin)<2
  error('must pass in more than one distribution')
end

groupLength=cellfun('length',varargin);

allData=cat(1,varargin{:});
ranks=crank(allData);

groupStartInAll=cumsum([1 groupLength]);

for groupNum=1:length(groupLength)
  groupIndices=groupStartInAll(groupNum)+(0:groupLength(groupNum)-1);

  groupRanks=ranks(groupIndices);
  
  C(groupNum)=sum(cos(groupRanks));
  S(groupNum)=sum(sin(groupRanks));  
end

W=2*sum((C.^2+S.^2)./groupLength);

p=1-chi2cdf(W,2*length(groupLength)-2); 
