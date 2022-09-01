function edf_split(eegfile,subject,outputdir)
% edf_split - Splits an EDF datafile into separate channels into
% the specified directory. 
%
% Acknowledgements: based on sopen.m in the biosig toolbox
%
% FUNCTION:
%    edf_split(eegfile,subject,outputdir)
% 
% INPUT ARGs:
% eegfile = 'UP014_23Oct07_1612.edf'
% subject = 'UP014'
% outputdir = '/data/eeg/UP014/eeg.noreref'
%

outdataformat = 'int16';

% check whether the output directory exists, and if not, create it
if ~exist(outputdir,'dir')
  mkdir(outputdir);
end

T0 = repmat(nan,1,6);

%%%%% Define Valid Data types %%%%%%
%GDFTYPES=[0 1 2 3 4 5 6 7 16 17 255+(1:64) 511+(1:64)];
GDFTYPES=[0 1 2 3 4 5 6 7 16 17 18 255+[1 12 22 24] 511+[1 12 22 24]];

%%%%% Define Size for each data type %%%%%
GDFTYP_BYTE=zeros(1,512+64);
GDFTYP_BYTE(256+(1:64))=(1:64)/8;
GDFTYP_BYTE(512+(1:64))=(1:64)/8;
GDFTYP_BYTE(1:19)=[1 1 1 2 2 4 4 8 8 4 8 0 0 0 0 0 4 8 16]';

H2idx = [16 80 8 8 8 8 8 80 8 32];        
ErrNo = 0; 
HANDEDNESS = {'unknown','right','left','equal'}; 
GENDER  = {'unknown','male','female'};
SCALE13 = {'unknown','no','yes'};
SCALE14 = {'unknown','no','yes','corrected'};


[fid,message] = fopen(eegfile,'rb','ieee-le');
if fid<0 
  ErrNo = [32,ErrNo];
  return;
end;

%%% read fixed header
[H1,count]=fread(fid,[1,192],'uchar');
if count<129
  ErrNo = [1,ErrNo];
  return;
end

versionNo=char(H1(1:8)); % 8 Byte  version number 
if ~(strcmp(versionNo,'0       '))
  ErrNo = [1,ErrNo];
  if ~strcmp(versionNo(1:3),'   '); % if not a scoring file, 
  %	    return; 
  end;
end;
if strcmp(char(H1(1:8)),'0       ') 
  versionNo= 0; 
elseif all(abs(H1(1:8))==[255,abs('BIOSEMI')]), 
  versionNo = -1; 
elseif all(H1(1:3)==abs('GDF'))
  versionNo = str2double(char(H1(4:8))); 
else
  ErrNo = [1,ErrNo];
  if ~strcmp(versionNo(1:3),'   '); % if not a scoring file, 
                                %	    return; 
  end;
end;

H1(193:256)= fread(fid,[1,256-192],'uchar');     %
H1 = char(H1); 
PID = deblank(char(H1(9:88)));  % 80 Byte local patient identification
RID = deblank(char(H1(89:168))); % 80 Byte local recording identification
[Patient.id,tmp] = strtok(PID,' ');
[tmp1,tmp] = strtok(tmp,' ');
[tmp1,tmp] = strtok(tmp,' ');
Patient.Name = tmp(2:end); 
recordDate = H1(169:176);
recordTime = H1(177:181);

tmp=(find((H1<32) | (H1>126))); 		%%% syntax for Matlab
if ~isempty(tmp) %%%%% not EDF because filled out with ASCII(0) - should be spaces
  ErrNo=[1025,ErrNo];
end;
                        
tmp = repmat(' ',1,22);
tmp([3:4,6:7,9:10,12:13,15:16,18:19]) = H1(168+[7:8,4:5,1:2,9:10,12:13,15:16]);
tmp1 = str2double(tmp);
if length(tmp1)==6,
  T0(1:6) = tmp1;
end;
                        
if any(isnan(T0)),
  ErrNo = [1032,ErrNo];
                                
  tmp = H1(168 + [1:16]);
  tmp(tmp=='.' | tmp==':' | tmp=='/' | tmp=='-') = ' ';
  tmp1 = str2num(tmp(1:8));
  if length(tmp1)==3,
    T0 = tmp1([3,2,1]);
  end;	
  tmp1 = str2num(tmp(9:16));
  if length(tmp1)==3,
    T0(4:6) = tmp1; 
  end;
  if any(isnan(T0)),
    ErrNo = [2,ErrNo];
  end;
end;
                        
% Y2K compatibility until year 2084
if T0(1) < 85    % for biomedical data recorded in the 1950's and converted to EDF
  T0(1) = 2000+T0(1);
elseif T0(1) < 100
  T0(1) = 1900+T0(1);
  %else % already corrected, do not change
end;
                        
HeadLen = str2double(H1(185:192));           % 8 Bytes  Length of Header
reserved1=H1(193:236);              % 44 Bytes reserved   
NRec    = str2double(H1(237:244));     % 8 Bytes  # of data records
Dur     = str2double(H1(245:252));     % 8 Bytes  # duration of data record in sec
NS      = str2double(H1(253:256));     % 4 Bytes  # of signals
AS.H1 = H1;	                     % for debugging the EDF Header
                        

if any(size(NS)~=1) %%%%% not EDF because filled out with ASCII(0) - should be spaces
  fprintf('Warning SOPEN (GDF/EDF/BDF): invalid NS-value in header of %s\n',eegfile);
  ErrNo=[1040,ErrNo];
  NS=1;
end;
% Octave assumes HDR.NS is a matrix instead of a scalare. Therefore, we need
% Otherwise, eye(HDR.NS) will be executed as eye(size(HDR.NS)).
NS = NS(1);     
                
if isempty(HeadLen) %%%%% not EDF because filled out with ASCII(0) - should be spaces
  ErrNo=[1056,ErrNo];
  HeadLen=256*(1+NS);
end;

if isempty(NRec) %%%%% not EDF because filled out with ASCII(0) - should be spaces
  ErrNo=[1027,ErrNo];
  NRec = -1;
end;
                
if isempty(Dur) %%%%% not EDF because filled out with ASCII(0) - should be spaces
  ErrNo=[1088,ErrNo];
  Dur=30;
end;
                
if  any(T0>[2084 12 31 24 59 59]) | any(T0<[1985 1 1 0 0 0])
  ErrNo = [4, ErrNo];
end;

idx1=cumsum([0 H2idx]);
idx2=NS*idx1;
                        
h2=zeros(NS,256);
[H2,count]=fread(fid,NS*256,'uchar');
if count < NS*256 
  ErrNo=[8,ErrNo];
  return; 
end;
                        
tmp = find((H2<32) | ((H2>126) & (H2~=255) & (H2~=181)& (H2~=230))); 
if ~isempty(tmp) %%%%% not EDF because filled out with ASCII(0) - should be spaces
  H2(tmp) = 32; 
  ErrNo = [1026,ErrNo];
end;

for k=1:length(H2idx);
  h2(:,idx1(k)+1:idx1(k+1))=reshape(H2(idx2(k)+1:idx2(k+1)),H2idx(k),NS)';
end;
h2=char(h2);

Label      =            h2(:,idx1(1)+1:idx1(2));
Transducer =            h2(:,idx1(2)+1:idx1(3));
PhysDim    =            h2(:,idx1(3)+1:idx1(4));
PhysMin    = str2num(h2(:,idx1(4)+1:idx1(5)));
PhysMax    = str2num(h2(:,idx1(5)+1:idx1(6)));
DigMin     = str2num(h2(:,idx1(6)+1:idx1(7)));
DigMax     = str2num(h2(:,idx1(7)+1:idx1(8)));
PreFilt    =            h2(:,idx1(8)+1:idx1(9));
AS.SPR     = str2num(h2(:,idx1(9)+1:idx1(10)));

% remove annotations from the Label and adjust the other variables
LabelC = cellstr(Label);
annotat1 = find(ismember(LabelC,'Pulse Rate'));
annotat2 = find(ismember(LabelC,'RR'));
annotat3 = find(ismember(LabelC,'IBI'));
annotat4 = find(ismember(LabelC,'Bursts'));
annotat5 = find(ismember(LabelC,'Suppr'));
hasAnnot = 0;
extraDat = 0; %if there are no annotation channels, there is no
              %extra data in each second of data samples
restInd = [];
if ~isempty(annotat1)
  restInd = [restInd annotat1];
end;
if ~isempty(annotat2)
  restInd = [restInd annotat2];
end
if ~isempty(annotat3)
  restInd = [restInd annotat3];
end;
if ~isempty(annotat4)
  restInd = [restInd annotat4];
end
if ~isempty(annotat5)
  restInd = [restInd annotat5];
end;
if ~isempty(restInd)
  restIndStart = min(restInd);
  % find out how many datapoints the extra annotation channels take
  % up (it is in position 217 in the description string)
  extraDat = sum(str2num(h2(restInd,217)));
  Label = Label(1:(restIndStart-1),:);
  NS = size(Label,1);
  chan = 1:NS;
  Transducer = Transducer(1:(restIndStart-1),:);
  PhysDim = PhysDim(1:(restIndStart-1),:);
  PhysMin = PhysMin(1:(restIndStart-1));
  PhysMax = PhysMax(1:(restIndStart-1));
  DigMin = DigMin(1:(restIndStart-1));
  DigMax = DigMax(1:(restIndStart-1));
  PreFilt = PreFilt(1:(restIndStart-1));
  AS.SPR = AS.SPR(1:(restIndStart-1));
  hasAnnot = 1;
end

if (versionNo ~= -1),
  GDFTYP     = 3*ones(1,NS);	%	datatype
else
  GDFTYP     = (255+24)*ones(1,NS);	%	datatype
end;

if isempty(AS.SPR), 
  fprintf('Warning SOPEN (GDF/EDF/BDF): invalid SPR-value in header of %s\n',eegfile);
  AS.SPR=ones(NS,1);
  ErrNo=[1028,ErrNo];
end;

%   threshold  = [DigMin,DigMax];       % automated overflow detection 
%   if (versionNo == 0) & FLAG.OVERFLOWDETECTION,   % in case of EDF and OVERFLOWDETECTION
%     fprintf('WARNING SOPEN(EDF): Physical Max/Min values of EDF data are not necessarily defining the dynamic range.\n'); 
%     fprintf('   Hence, OVERFLOWDETECTION might not work correctly. See also EEG2HIST and read \n'); 
%     fprintf('   http://dx.doi.org/10.1016/S1388-2457(99)00172-8 (A. Schlˆgl et al. Quality Control ... Clin. Neurophysiol. 1999, Dec; 110(12): 2165 - 2170).\n'); 
%     fprintf('   A copy is available here, too: http://www.dpmi.tugraz.at/schloegl/publications/neurophys1999_2165.pdf \n'); 
%   end;
% end; 
  
if any(PhysMax==PhysMin), ErrNo=[1029,ErrNo]; end;	
if any(DigMax ==DigMin ), ErrNo=[1030,ErrNo]; end;	
Cal = (PhysMax-PhysMin)./(DigMax-DigMin);
Off = PhysMin - Cal .* DigMin;
  
AS.SampleRate = AS.SPR / Dur;
SPR=1;
chan = 1:NS;
if strcmp(reserved1(1:4),'EDF+')
  tmp = strmatch('EDF Annotations',Label);
  chan(tmp)=[];
end; 
  
for k=chan,
  if (AS.SPR(k)>0)
    SPR = lcm(SPR,AS.SPR(k));
  end;
end;
SampleRate = SPR/Dur;
  
AS.spb = sum(AS.SPR);	% Samples per Block
AS.bi  = [0;cumsum(AS.SPR(:))]; 
AS.BPR = ceil(AS.SPR.*GDFTYP_BYTE(GDFTYP+1)'); 
while any(AS.BPR ~= AS.SPR.*GDFTYP_BYTE(GDFTYP+1)');
 fprintf('\nError SOPEN (GDF/EDF/BDF): block configuration in file %s not supported.\n',eegfile);
end;
AS.SAMECHANTYP = all(AS.BPR == (AS.SPR.*GDFTYP_BYTE(GDFTYP+1)')) & ~any(diff(GDFTYP)); 
AS.bpb = sum(ceil(AS.SPR.*GDFTYP_BYTE(GDFTYP+1)'));	% Bytes per Block
Calib  = [Off'; diag(Cal)];


% filesize, position of eventtable, headerlength, etc. 	
AS.EVENTTABLEPOS = -1;
w = dir(eegfile);
AS.endpos = w.bytes;
if (AS.endpos == HeadLen)
  NRec = 0; 
elseif NRec == -1   % unknown record size, determine correct NRec
  NRec = floor((AS.endpos - HeadLen) / AS.bpb);
end
if  (NRec*AS.bpb) ~= (AS.endpos - HeadLen);
  ErrNo= [16,ErrNo];
  tmp = NRec; 
  NRec = floor((AS.endpos - HeadLen) / AS.bpb);
  if tmp~=NRec,
    fprintf('\nWarning SOPEN (EDF/BDF): filesize (%i) of %s does not fit headerinformation (NRec = %i not %i).\n',AS.endpos,eegfile,tmp,NRec);
    NRec = tmp; % adapt NRec to reflect actual # of samples
  else
    fprintf('\nWarning: incomplete data block appended (ignored) in file %s.\n',eegfile);
  end
else
  AS.EVENTTABLEPOS = HeadLen + AS.bpb*NRec;
end;

% prepare SREAD for different data types 
n = 0; 
typ = [-1;GDFTYP(:)];
for k = 1:NS; 
  if (typ(k) == typ(k+1)),
    AS.c(n)   = AS.c(n)  + AS.SPR(k);
    AS.c2(n)  = AS.c2(n) + AS.BPR(k);
  else
    n = n + 1; 
    AS.c(n)   = AS.SPR(k);
    AS.c2(n)  = AS.BPR(k);
    AS.TYP(n) = GDFTYP(k);
  end;
end;


%%% write params.txt file & jacksheet
allMonths = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
monthStr = allMonths{str2num(recordDate(4:5))};
filestem = fullfile(outputdir,[subject '_' recordDate(1:2) monthStr recordDate(7:8) '_' recordTime(1:2) recordTime(4:5)]);
[pathstr,name,ext] = fileparts(filestem);
paramfile = fullfile(pathstr,'params.txt');
fout1 = fopen(paramfile,'w','l');
%% FIX THE GAIN!!!
fprintf(fout1,'samplerate %.2f\ndataformat ''%s''\ngain %g\n',SPR,outdataformat,Cal(1));
fclose(fout1);
% new params matching base name
paramfile = [filestem '.params.txt'];
fout2 = fopen(paramfile,'w');
fprintf(fout2,'samplerate %d\ndataformat ''%s''\ngain %g\n',SPR,outdataformat,Cal(1));
fclose(fout2);

%print the labels for each channel
fout3 = fopen(fullfile(outputdir,'jacksheet.txt'),'w','l');
for c = 1:size(Label,1)
  fprintf(fout3,'%d %s\n',c,Label(c,:));
end
fclose(fout3);

% start reading the data
% (it is specified as:
% nr of samples[1] * integer : first signal in the data record 
% nr of samples[2] * integer : second signal 
% .. 
% .. 
% nr of samples[ns] * integer : last signal

% go to the correct position in the file
status = fseek(fid, HeadLen, 'bof');
FILE.POS  = 0;
FILE.OPEN = 1;

% read the data from the file
fprintf('reading samples...\n');
datLen = NRec*SPR;
raw = int16(zeros(NS,datLen));
for t = 1:NRec
  % the factor extraDat deals with the annotations
  dat = fread(fid,NS*SPR + extraDat,'int16');
  dat = dat(1:(end-extraDat));
  if length(dat)>0
    if length(dat)<NS*SPR
      fprintf('too long\n');
      keyboard
    else
      dat2= reshape(dat,SPR,NS);
      raw(:,(t-1)*SPR+1:t*SPR) = int16(dat2');
    end
  end
end
  
% loop over channels
fprintf('Processing %d channels:\n',length(chan));
for c = chan
  fprintf('%d ',c);
  %%% write channel
  % make the filename
  chanfile = sprintf('%s.%03i', filestem,c);  
  % open and write the file
  fchan = fopen(chanfile,'w','l');
  fwrite(fchan,raw(c,:),outdataformat);
  fclose(fchan);
end
fprintf('\n');
