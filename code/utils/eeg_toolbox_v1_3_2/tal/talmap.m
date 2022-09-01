function [positives,negatives]=talmap(maptype,layoutnum,subcol,subshape,funfstem,pthresh,excludecriteria)
% [positives,negatives]=talmap(maptype,layoutnum,subcol,subshape,funfstem,pthresh,excludecriteria)
%
% maptype: 'num' - lead # labels [put [] for subshape & funfstem & pthresh]
%          'mrkfilled' - plot a filled shape marker for each lead
%                        [put [] for subshape & funfstem & pthresh]
%          'mrkempty'  - plot an empty shape marker for each lead
%                        [put [] for subshape & funfstem & pthresh]
%          '-+thresh'   - A thresholded map for p-value thresholds:
%                        |value| > pthresh - not plotted
%                        0 < value < pthresh   - filled subcol
%                        -pthresh < value < 0  - filled with inverse of subcol
%                      This is the most useful map type
%          '+-thresh'   - A thresholded map:
%                        |value| < pthresh - not plotted
%                        value > pthresh  - filled subcol
%                        value < -pthresh - filled with inverse of subcol
%                      This is the most useful map type
%          'fun'     - functional map, intensity reflects value.
%                      Uses the inverse colour for negative values
%                      *** pthresh is a 2-vector: [minval maxval]
%          '-logp'   - functional map using the sign(p)*(-log(p)) transform
%***          'phase'   - a phase symbol map
%***          'plf'     - a phase-locking factor clock-hand map
%***          'fourier' - a Fourier component clock-hand map
% excludecriteria: 0 - don't exclude; plot everything
%                  1 - exclude; only plot included electrodes
%                  2 - only plot excluded electrodes
%
% Requires a complete tal/tal_params.txt file
% retrieve the positions of the views

okay=0;

if(strcmp(maptype,'num')) okay=1; end;
if(strcmp(maptype,'mrkfilled')) okay=1; end;
if(strcmp(maptype,'mrkempty')) okay=1; end;
if(strcmp(maptype,'-+thresh')) okay=1; end;
if(strcmp(maptype,'+-thresh')) okay=1; end;
if(strcmp(maptype,'fun')) 
  if(length(pthresh)~=2)
    fprintf(1,'\n\nFatal Error: Type ''fun'' needs a 2-vector for pthresh: pthresh=[minval maxval]\n\n');
    return
   %break;
  else okay=1; end;
end;
if(strcmp(maptype,'-logp')) okay=1; end;
if(strcmp(maptype,'phase')) okay=1; end;
if(strcmp(maptype,'plf')) okay=1; end;
if(strcmp(maptype,'fourier')) okay=1; end;

if(~okay)
   fprintf(1,'\nFatal Error: I don''t know maptype "%s"\n\n',maptype);
   return
   %break;
end

xc=-999;

[views,viewpos,origins,scales,xcol,ycol,xsign,ysign,XY]=get_tal_params(layoutnum);

[whichviews,viewpos]=tal_view(0,layoutnum);

new_coords_inf=get_tal_montage('inf');
new_coords_lsag=get_tal_montage('lsag');
new_coords_rsag=get_tal_montage('rsag');
new_coords_lsag_int=get_tal_montage('lsag_int');
new_coords_rsag_int=get_tal_montage('rsag_int');
new_coords_hipp=get_tal_montage('hipp');

ncoords={new_coords_inf,new_coords_lsag,new_coords_rsag,new_coords_lsag_int,new_coords_rsag_int,new_coords_hipp};



if(~exist('excludecriteria')) excludecriteria=1; end;

% Settings that rely on eeganalparams

forgrant=eeganalparams('forgrant'); forslides=eeganalparams('forslides');
if(~exist('lead_rad')) lead_rad=eeganalparams('lead_rad'); end;
if(~exist('excludeleads')) excludeleads=eeganalparams('excludeleads'); end;
edgethresh=0.9; % threshold brightness to put a border on ('hotfun') only
if(forslides & istemplate) edgec=[1 1 1]; else edgec=[]; end;
excludeval=999;
vectorscale=4; % scale factor for drawing phase/fourier/plf vectors

% ensure that the brightness is the same
brightness=norm(subcol); subcolneg=[1 1 1]-subcol;
subcolpos=subcol; if(subcol==[0 0 0]) subcolneg=subcol;
else subcolneg=0.75*subcolneg*brightness/norm(subcol); end

% three types of maps: numerical, shape, circle or vector
if(strcmp(maptype,'num')) mtypeno=0;
elseif( strcmp(maptype,'phase') | strcmp(maptype,'plf') | strcmp(maptype,'fourier') ) mtypeno=2;
else mtypeno=1; end;

% perform the exclude/include criteria
switch(excludecriteria)
  case 0
  case 1
    goodleads=getleads('good_leads.txt');
    for V=whichviews
      if(~isempty(ncoords{V}))
        [values,keep]=intersect(ncoords{V}(:,1),goodleads);
	ncoords{V}=ncoords{V}(keep,:);
      end
    end % view loop, excluding bad leads
  case 2
    goodleads=getleads('good_leads.txt');
    for V=whichviews
      if(~isempty(ncoords{V}))
        [values,keep]=setdiff(ncoords{V}(:,1),goodleads);
	ncoords{V}=ncoords{V}(keep,:);
      end
    end % view loop, excluding _good_ leads
end

positives=[]; negative=[];

% Create the electrode values, for subsequent plotting

istemplate=0;

leads=getleads('leads.txt');

if(strcmp(maptype,'num')) % nothing left to do
elseif(strcmp(maptype,'mrkfilled'))
  func=[leads;ones(1,length(leads))]; allcol=ones(length(leads),1)*subcol;
elseif(strcmp(maptype,'mrkempty')) % use white colour for all
  func=[leads;ones(1,length(leads))]; allcol=ones(length(leads),1)*[1 1 1];
  istemplate=1;
else % read in the functional file
  if(~isstr(funfstem)) func=funfstem;
  else
    funfname=sprintf('fun/%s.fun',funfstem);
    if(~exist(funfname))
      fprintf(1,'\nError: Can''t open functional file ''%s''\n\n',funfname);
      return
      %break;
    end
    in=fopen(funfname,'r');
    if(mtypeno==2)
      func=fscanf(in,'%i %f %f',[3,inf]);
      func=[func(1,:);func(2,:)+i*func(3,:)];
    else func=fscanf(in,'%i %f',[2,inf]); end;   fclose(in);
  end
end % read in the functional file

if(strcmp(maptype,'fun'))
  positives=find(func(2,:)>=0);
  if(~isempty(positives)) 
    minfun=pthresh(1); maxfun=pthresh(2); rangefun=maxfun-minfun;
    thecolour=(func(2,positives)-minfun)/rangefun;
    thecolour(find(thecolour>1))=1; thecolour(find(thecolour<0))=0;
    allcol(positives,:)=thecolour'*subcolpos;
  end
  negatives=find(func(2,:)<0); 
  if(~isempty(negatives)) 
    thecolour=(abs(func(2,negatives))-minfun)/rangefun;
    thecolour(find(thecolour>1))=1; thecolour(find(thecolour<0))=0;
    allcol(negatives,:)=thecolour'*subcolneg;
  end
elseif(strcmp(maptype,'-logp'))
  minP=-log10(0.05);
  positives=find(func(2,:)>=0); negatives=find(func(2,:)<0);
  thecolour=-log10(abs(func(2,:))); thecolour=(thecolour-minP)/(pthresh-minP);
  thecolour(find(thecolour>1))=1; thecolour(find(thecolour<0))=0;
  if(~isempty(positives))
    allcol(positives,:)=thecolour(positives)'*subcolpos;
  end
  if(~isempty(negatives))
    allcol(negatives,:)=thecolour(negatives)'*subcolneg;
  end
elseif(strcmp(maptype,'-+thresh'))
  positives=find( (func(2,:)<=pthresh) & (func(2,:)>=0) );
  negatives=find( (func(2,:)>=(-pthresh)) & (func(2,:)<0) );
  unpasses=setdiff(setdiff(1:size(func,2),positives),negatives);
  if(~isempty(positives))
    allcol(positives,:)=ones(size(positives))'*subcolpos;
  end
  if(~isempty(negatives))
    allcol(negatives,:)=ones(size(negatives))'*subcolneg;
  end
  if(~isempty(unpasses)) allcol(unpasses,:)=xc*ones(length(unpasses),3);end;
elseif(strcmp(maptype,'+-thresh'))
  positives=find(func(2,:)>=pthresh); negatives=find(func(2,:)<=(-pthresh));
  unpasses=setdiff(setdiff(1:size(func,2),positives),negatives);
  if(~isempty(positives))
    allcol(positives,:)=ones(size(positives))'*subcolpos;
  end
  if(~isempty(negatives))
    allcol(negatives,:)=ones(size(negatives))'*subcolneg;
  end
  if(~isempty(unpasses)) allcol(unpasses,:)=xc*ones(length(unpasses),3); end
end % assigning values based on maptype


% ***** Finally, plot the leads on each of the views *****
for V=whichviews
  if(~isempty(ncoords{V}))
    nc=ncoords{V}; sc=scales(V); O=origins{V}+viewpos{V};
    xi=xcol(V); yi=ycol(V); xfac=xsign(V)*sc; yfac=ysign(V)*sc;
    for l=1:size(nc,1)
      XXX=O(1)+xfac*nc(l,xi); YYY=O(2)+yfac*nc(l,yi);
      if(mtypeno==0)
        h=text(XXX,YYY,sprintf('%i',nc(l,1)));
        set(h,'FontWeight','Bold','Color',subcol,'HorizontalAlignment','center');
      elseif(strcmp(maptype,'phase')) % phase map
        thelead=find(func(1,:)==nc(l,1));
        if(~isempty(thelead))
          drawshape('phase',XXX,YYY,vectorscale*lead_rad*sc,subcol,func(2,thelead),edgec);
        end
      elseif(mtypeno==2)
        thelead=find(func(1,:)==nc(l,1));
        if(~isempty(thelead))
	  drawshape('phase',XXX,YYY, ...
              norm(func(2,thelead))*vectorscale*lead_rad*sc,subcol, ...
              angle(func(2,thelead)),edgec);
        end;
        thelead=find(func(1,:)==nc(l,1));
        if(~isempty(thelead))
	  if(exist('allcol'))
            if((norm(allcol(thelead,:))/norm([1 1 1]))>edgethresh) EDGEC=[0 0 0];
	    else EDGEC=edgec; end;
	      if(istemplate)
	         drawshape(subshape,XXX,YYY, ...
                lead_rad*sc,[0 0 0],edgec);
	      elseif(allcol(thelead,1)~=xc)
		drawshape(subshape,XXX,YYY,lead_rad*sc,allcol(thelead,:),EDGEC);
	      end
            end
         end % if allcol exists
      else % mtypeno==1
        thelead=find(nc(l,1)==func(1,:));
        if(~isempty(thelead))
          if((norm(allcol(thelead,:))/norm([1 1 1]))>edgethresh) EDGEC=[0 0 0];
          else EDGEC=edgec; end;
          if(istemplate)
            drawshape(subshape,XXX,YYY,lead_rad*sc,[0 0 0],edgec);
          elseif(allcol(thelead,1)~=xc)
            drawshape(subshape,XXX,YYY,lead_rad*sc,allcol(thelead,:),edgec);
          end
	end % if there are leads
      end % all the functional maps
    end % lead loop
  end % if there are leads on this view
end % looping through views

if(exist('positives') & exist('func') )
  positives=func(1,positives); positives=intersect(positives,goodleads);
end;
if(exist('negatives') & exist('func') )
  negatives=func(1,negatives); negatives=intersect(negatives,goodleads);
end;
