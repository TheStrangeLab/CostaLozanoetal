function modes=jemd(data)
numModes=10;
numSifts=10;

modes=[];

x=1:length(data);
      

residue=data;
n=10; plotNum=2;

while ~doneSifting(residue)
  [imf,residue]=doSift(residue);
  modes=[modes;imf];
  
  plotNum=plotNum+1;
end

modes=[modes;residue];

function y=doneSifting(d)
%We are done sifting is there a monotonic function
y=sum(localmax(d))+sum(localmax(-d))<=2;


function [imf,residue]=doSift(data)
%This function is modified to use the sifting-stopping criteria from Huang et al
%(2003) (this is the suggestion of Peel et al., 2005).  Briefly, we sift until the
%number of extrema and zerocrossings differ by at most one, then we continue
%sifting until the number of extrema and ZCs both remain constant for at least
%five sifts.

imf=data;

while 1 %sift until num extrema and ZC differ by at most 1
  imf=doOneSift(imf);
  [numExtrema,numZC]=analyzeIMF(imf);
  %disp(sprintf('numextrema=%d, numZC=%d', numExtrema, numZC)); 
  if abs(numExtrema-numZC)<=1
    break;
  end
end

%then continue until numExtrema and ZCs are constant for at least 5 sifts
%(Huang et al., 2003)

numConstant=5; 
while 1
  imf=doOneSift(imf);
  [numExtrema(end+1),numZC(end+1)]=analyzeIMF(imf);
  %disp(sprintf('FINAL STAGE: numextrema=%d, numZC=%d', numExtrema(end), numZC(end))); 
  if length(numExtrema)>=numConstant && ...
        all(numExtrema(end-4:end)==numExtrema(end)) && ...
        all(numZC(end-4:end)==numZC(end))
    break
  end
end


residue=data-imf;


function detail=doOneSift(data)

upper=getUpperSpline(data);
lower=-getUpperSpline(-data);
%upper=jinterp(find(maxes),data(maxes),xs);
%lower=jinterp(find(mins),data(mins),xs);

imf=mean([upper;lower],1);
detail=data-imf;

%plot(xs,data,'b-',xs,upper,'r--',xs,lower,'r--',xs,imf,'k-')


function s=getUpperSpline(data)
maxInds=find(localmax(data));

%this is the Mirroring algoirthm from Rilling et al. (2003)

if length(maxInds)==1
  %Special case: if there is just one max, then entire spline is that number  
  s=repmat(data(maxInds),size(data));
  return
end

%Start points
if maxInds(1)==1   %if first point is a local max
  preTimes=2-maxInds(2);
  preData=data(maxInds(2));
else %if first point is NOT local max
  preTimes=2-maxInds([2 1]);
  preData=data(maxInds([2 1]));
end

%end points
if maxInds(end)==length(data) %if last point is a local max
  postTimes=2*length(data)-maxInds(end-1);
  postData=data(maxInds(end-1));
else %if last point is NOT a local max
  postTimes=2*length(data)-maxInds([end end-1]);
  postData=data(maxInds([end end-1]));
end

t=[preTimes maxInds postTimes];
d2=[preData data(maxInds) postData];
s=interp1(t,d2,1:length(data),'spline');

%plot(1:length(data),data,'b-',1:length(data),s,'k-',t,d2,'r--');
%keyboard
  

function [numExtrema,numZC]=analyzeIMF(d)
xs=1:length(d);

numExtrema=sum(localmax(d))+sum(localmax(-d));
numZC=sum(diff(sign(d))~=0);

% if debug
%   clf
%   a1=subplot(2,1,1);
%   plot(xs,d,'b-',xs,upper,'k-',xs,lower,'k-');
%   axis tight;
  
%   a2=subplot(2,1,2);
%   plot(xs,stopScore,'b-',[0 length(d)],[thresh1 thresh1],'k--',[0 length(d)],[thresh2 ...
%                       thresh2],'r--');
%   axis tight;
%   xlabel(sprintf('score = %.3g',s));  
%   linkaxes([a1 a2],'x')
%   keyboard
  
% end






function yi=jinterp(x,y,xi);
if length(x)==1
  yi=repmat(y,size(xi));
else
  yi=interp1(x,y,xi,'spline');
end

  



function maxima=localmax(d)
%d must be a row vector

diffScore=diff(sign(diff([-inf d -inf])));
%this gets a value of -2 if it is an unambiguous local max
%value -1 denotes that the run its a part of may contain a local max


%Run length code with help from:
%  http://home.online.no/~pjacklam/matlab/doc/mtt/index.html
%(this is all painfully complicated, but I did it in order to avoid loops...)

%here calculate the position and length of each run
runEndingPositions=[find(d(1:end-1)~=d(2:end)) length(d)];
runLengths=diff([0 runEndingPositions]); %length of each consecutive run
runStarts=runEndingPositions-runLengths+1;

%Now concentrate on only the runs with length>1
realRunStarts=runStarts(runLengths>1);
realRunStops=runEndingPositions(runLengths>1);
realRunLengths=runLengths(runLengths>1);

%save only the runs that are local maxima
maxRuns=diffScore(realRunStarts)==-1 & diffScore(realRunStops)==-1;

%If a run is a local max, then count the middle position (rounded) as the 'max'
maxRunMiddles=round(realRunStarts(maxRuns)+realRunLengths(maxRuns)/2-1);

maxima=diffScore==-2;
maxima(maxRunMiddles)=true;

%make sure beginning & end are not local maxes
%maxima([1 end])=false;
  



