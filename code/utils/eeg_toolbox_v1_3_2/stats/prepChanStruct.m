% prepare a structure with channels for use with fieldTrip
function [chanStr,labels] = prepChanStruct

chanFile = 'GSN129c.sfp';

[chanID,x,y,z] = textread(chanFile,'%s\t%f\t%f\t%f');

chanStr.label = cell(length(x),1);
chanStr.pnt = zeros(length(x),3);
labels =cell(length(x),1);

% fill the channel struct
for e = 1:length(x)
  chanStr.label{e} = chanID(e);
  labels(e) = chanID(e);
  chanStr.pnt(e,:) = [x(e) y(e) z(e)];
end

allChans = 1:length(x);
badChans= [127 126 17 128 125 120 44 49 56 63 69 74 82 89 95 100 108 114];
goodChans = setdiff(allChans,badChans);

labels = labels(goodChans);
chanStr.label = chanStr.label(goodChans);
chanStr.pnt = chanStr.pnt(goodChans,:);
