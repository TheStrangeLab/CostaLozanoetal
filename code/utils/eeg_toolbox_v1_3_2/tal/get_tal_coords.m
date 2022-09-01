function new_coords=get_tal_coords(talfile,selectthese)
% new_coords=get_tal_coords(talfile,[selectthese])
%
% new_coords- the Talairach coordinates        - lead# inf-sup ant-post L-R
%
% This function converts raw coordinates into Talairach coordinates.
%
% To read, do:
% x=get_tal_coords('raw_coords.txt');
% fprintf(1,'%i: (%.1f,%.1f,%.1f)\n',x')
%
% To specify which electrode coordinates you want, optionally pass a
% vector containing those electrode numbers in the second parameter,
% selectthese.
%
% NOTE: Use uppercase talfile name (e.g., 'LSAG_coords.txt' instead
% of 'lsag_coords.txt') if the coordinates in the talfile are absolute
% instead of relative to some coregistration
%
% In the tal_params.txt file, specify the Talairach registration type: 
%  '3D' - from 3-D coordinate space, where you have electrodes and
%         landmarks all in the same 3D coord. system.
%  'slice' - from lots of slices, where electrodes are all on axial
%            slices, and you coregister and use sagittal and axial MRs to
%            obtain Talairach landmarks. This was the old method.



ind = 1;
if ~exist(talfile,'file')
  talfile = fullfile('tal',talfile);
  ind = ind+4;
end

if(talfile(ind)<'a')&exist(talfile,'file') 
  % no need to convert, so just set dummy TP
  TP.taltype = '3D';
  absolute=1; 
else 
  % will convert, so get all the parameters
  TP=get_tal_par;

  absolute=0; 
end

if ~absolute
  % The standardised Talairach constants
  Tal_inf_sup=67; % mm
  Tal_AC_inf_sup=80; % mm. This is the AC "origin"
  inf_sup_ratio=Tal_inf_sup/norm(TP.AC-TP.VAC_sup);
  Tal_ant=64; % mm
  Tal_post=96; % mm
  ant_ratio=Tal_ant/norm(TP.AC-TP.ACPC_ant);
  post_ratio=Tal_post/norm(TP.AC-TP.ACPC_post);
  Tal_axial_half=58; % mm
end

if(strcmp(TP.taltype,'3D'))
  in=fopen(sprintf('%s',talfile),'r');
  if( (in==-1) & (absolute) )
    absolute=0; in=fopen(sprintf('%s',lower(talfile)),'r');
  end
  if(in==-1)
    fprintf(1,'\nError: I couldn''t open the file %s\n',talfile);
    new_coords=[];
  elseif(absolute) % upper-case; read in absolute Talairach coords.
    new_coords=fscanf(in,'%f',[4,inf])'; fclose(in); % read coords directly
						     % columns:    [electrode #s]  [R-L]  [ant-post]  [inf-sup]
  else % relative to a coregistration
    raw_coords=fscanf(in,'%f',[4,inf]); fclose(in); % read in raw coords.
    new_coords(:,1)=raw_coords(1,:)'; nLeads=size(new_coords,1);
    
    % Set up the projections
    acpc=TP.PC-TP.AC;
    % A-P
    % Project ACPC_post onto ACPC
    proj=(TP.ACPC_post-TP.AC)*(acpc)'/(norm(acpc)^2); acpc_post=TP.AC + proj*acpc;
    proj=(TP.ACPC_ant-TP.AC)*(acpc)'/(norm(acpc)^2); acpc_ant=TP.AC + proj*acpc;
    
    for l=1:nLeads
      xyz=raw_coords(2:4,l)';
      % X (L-R)          Two cases
      proj=(xyz-TP.AC)*(TP.AC_left-TP.AC)'/(norm(TP.AC_left-TP.AC)^2);
      if(sign(proj)==1) % left side
	new_coords(l,2)=-proj;
      else % right side
	new_coords(l,2)=(xyz-TP.AC)*(TP.AC_right-TP.AC)'/...
	    (norm(TP.AC_right-TP.AC)^2);
      end
      new_coords(l,2)=new_coords(l,2)*Tal_axial_half;
      
      % y (A-P)          Three cases
      proj=(xyz-TP.AC)*(acpc_ant-TP.AC)'/(norm(acpc_ant-TP.AC)^2);
      if(sign(proj)==1) % anterior
	new_coords(l,3)=proj*Tal_ant;
      else
	proj=(xyz-TP.PC)*(-acpc)'/(norm(acpc)^2);
	if(sign(proj)==1) % along the AC-PC line
	  new_coords(l,3)=-((xyz-TP.AC)*acpc'/(norm(acpc)^2))*Tal_ant;
	else % posterior
	  new_coords(l,3)=-((xyz-TP.AC)*(acpc_post-TP.AC)'...
			   /(norm(acpc_post-TP.AC)^2))*Tal_post;
	end % posterior
      end % not ant
      
      % z (I-S)          One case
      new_coords(l,4)=(xyz-TP.PC)*(TP.VPC_sup-TP.PC)'/...
	  (norm(TP.VPC_sup-TP.PC)^2);
      new_coords(l,4)=new_coords(l,4)*Tal_inf_sup;
    end % looping through leads
  end % Can open the file; doing the transformation
else % taltype 'slice'
  
  % Are we going to try to be in absolute coords? The convention is
  % files containing absolute coordinates start with an uppercase
  % letter.
  %if(talfile(1)<'a') absolute=1; else absolute=0; end;
  
  % Read in the raw coordinate values.....
  % Col 1: lead #s; Col 2: slice; Col 3: x in-plane; Col 4: y in-plane
  % Col 5: x sagittal; Col 6: y sagittal
  in=fopen(sprintf('%s',talfile),'r');
  if( (in==-1) & (absolute) )
    absolute=0; in=fopen(sprintf('%s',lower(talfile)),'r');
  end
  if(in==-1)
    fprintf(1,'\nError: I couldn''t open the file %s\n',talfile);
    new_coords=[];
  elseif(absolute) % upper-case; read in absolute Talairach coords.
    new_coords=fscanf(in,'%f',[4,inf])'; fclose(in); % read coords directly
    % columns:    [electrode #s]  [R-L]  [ant-post]  [inf-sup]
  else % relative to a coregistration
      axial_L_ratio=...
      Tal_axial_half/mean([norm(TP.HAC_left-TP.HAC) norm(TP.HPC-TP.HPC_left)]);
  axial_R_ratio=...
      Tal_axial_half/mean([norm(TP.HAC_right-TP.HAC) norm(TP.HPC-TP.HPC_right)]);
  axial_half_ratio=mean([axial_L_ratio axial_R_ratio]);
  

    raw_coords=fscanf(in,'%f',[6,inf]); fclose(in); % read in raw coords.
    
    % lead numbers
    new_coords(:,1)=raw_coords(1,:)'; nLeads=size(new_coords,1);
    
    % axial coords: L-R (x-axis)
    for l=1:nLeads % do the slice coords. for each electrode
      sliceant(l,:)=TP.Slice_ant(find(TP.SlicesUsed==raw_coords(2,l)),:);
      slicepost(l,:)=TP.Slice_post(find(TP.SlicesUsed==raw_coords(2,l)),:);    
      slicenorm(l)=norm(sliceant(l,:)-slicepost(l,:));
      rawnorm(l)=norm(raw_coords(3:4,l)'-slicepost(l,:));
    end % lead loop
    midline=sliceant-slicepost; apangle=angle(midline(:,1)+i*midline(:,2));
    % Project onto the midline, then compute the dot product (adjacent side)
    y=sum((raw_coords(3:4,:)'-slicepost).*(midline),2); y=y./slicenorm';
    % Pythagoras: rawnorm is the hypotenuse, y is the adjacent side--
    % projection onto the midline
    z=sqrt(rawnorm'.^2 - y.^2);
    X=raw_coords(3:4,:)'-slicepost; Xangle=angle(X(:,1)+i*X(:,2));
    these=find(Xangle<apangle); % Right hemisphere
    new_coords(these,2)=z(these)*axial_R_ratio;
    these=find(Xangle>=apangle); % Left hemisphere
    new_coords(these,2)=-z(these)*axial_L_ratio;
    
    % sagittal coords: ant-post (y-axis)
    val=(raw_coords(5:6,:)'+ones(nLeads,1)*TP.ACPC_post)*(TP.ACPC_ant-TP.ACPC_post)'/norm(TP.ACPC_ant-TP.ACPC_post)*post_ratio;
    these=(find(val>Tal_post)); % in the posterior region
    new_coords(these,3)=0+(raw_coords(5:6,these)'-ones(length(these),1)*TP.AC)*(TP.ACPC_ant-TP.AC)'/norm(TP.ACPC_ant-TP.AC)*post_ratio;
    these=(find(val<=Tal_post)); % in the anterior region
    new_coords(these,3)=0+(raw_coords(5:6,these)'-ones(length(these),1)*TP.AC)*(TP.ACPC_ant-TP.AC)'/norm(TP.ACPC_ant-TP.AC)*ant_ratio;
    
    % sagittal coords: inf-sup (z-axis)
    new_coords(:,4)=(raw_coords(5:6,:)'-ones(nLeads,1)*TP.AC)*(TP.VAC_sup-TP.AC)'/norm(TP.VAC_sup-TP.AC)*inf_sup_ratio;
  end % if the talfile exists
  
end % 'slice' type 

% If we were passed a set of electrodes, then scrape off these
if(exist('selectthese'))
  new_coords=new_coords(find(ismember(new_coords(:,1),selectthese)),:);
end
