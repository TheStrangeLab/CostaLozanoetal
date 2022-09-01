function TP=get_tal_par;
% TP=get_tal_par;

if(nargin<1) layoutnum=0; end;

% The standardised Talairach constants
TP.Tal_inf_sup=75; % mm
TP.Tal_sup=TP.Tal_inf_sup; TP.Tal_inf=TP.Tal_inf_sup/(1.7);
TP.Tal_ant=70; % mm
TP.Tal_post=103; % mm
TP.Tal_axial_half=61; % mm

if exist('tal_params.txt','file')
  in = fopen('tal_params.txt','r');
elseif exist('tal/tal_params.txt','file')
  in = fopen('tal/tal_params.txt','r');
else
 fprintf(1,'\nError: couldn''t open file tal_params.txt\n');
 return;
end

taltype=fscanf(in,'%s',1);
if(strcmp(lower(taltype(1:2)),'3d'))
  % ***************** Do the 3D coordinates version *****************
  TP.taltype='3D';
  field=fscanf(in,'%s',1);
  while(~isempty(field))
    switch(lower(field))
     case 'ac', TP.AC=fscanf(in,'%f',[1,3]);
     case 'pc', TP.PC=fscanf(in,'%f',[1,3]);
     case 'acpc_post', TP.ACPC_post=fscanf(in,'%f',[1,3]);
     case 'acpc_ant', TP.ACPC_ant=fscanf(in,'%f',[1,3]);
     case 'vac_sup', TP.VAC_sup=fscanf(in,'%f',[1,3]);
     case 'vpc_sup', TP.VPC_sup=fscanf(in,'%f',[1,3]);
     case 'ac_right', TP.AC_right=fscanf(in,'%f',[1,3]);
     case 'ac_left', TP.AC_left=fscanf(in,'%f',[1,3]);
     case 'pc_right', TP.PC_right=fscanf(in,'%f',[1,3]);
     case 'pc_left', TP.PC_left=fscanf(in,'%f',[1,3]);
     case 'inf_trace_ac', TP.inf_trace_AC=fscanf(in,'%i',[1,2]);
     case 'inf_trace_pc', TP.inf_trace_PC=fscanf(in,'%i',[1,2]);
     case 'inf_trace_scale', TP.inf_trace_scale=fscanf(in,'%f',1);
     case 'l_trace_ac', TP.L_sag_trace_AC=fscanf(in,'%i',[1,2]);
     case 'l_trace_pc', TP.L_sag_trace_PC=fscanf(in,'%i',[1,2]);
     case 'r_trace_ac', TP.R_sag_trace_AC=fscanf(in,'%i',[1,2]);
     case 'r_trace_pc', TP.R_sag_trace_PC=fscanf(in,'%i',[1,2]);
     case 'sag_trace_scale', TP.sag_trace_scale=fscanf(in,'%f',1);
     otherwise
      fprintf(1,'\nError: I don''t recognise the field %s\n',field);
    end % switch
    field=fscanf(in,'%s',1);
  end % while reading through the file
  fclose(in);

  % Error messages
  fld='AC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='PC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='ACPC_post'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='ACPC_ant'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='VAC_sup'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='VPC_sup'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='AC_right'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='AC_left'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='PC_right'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='PC_left'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='inf_trace_AC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='inf_trace_PC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='inf_trace_scale'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='L_sag_trace_AC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='L_sag_trace_PC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='R_sag_trace_AC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='R_sag_trace_PC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='sag_trace_scale'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  

else % ******************* Do the Slice version *********************
  TP.taltype='slice';
  fseek(in,0,-1);
  %fclose(in); in=fopen('tal/tal_params.txt','r');
  TP.Slice_ant=[]; TP.Slice_post=[]; TP.SlicesUsed=[];
  field=fscanf(in,'%s',1);
  while(~isempty(field))
    if( (field(1)>='0') & (field(1)<='9')) % if a number, do a slice
      snum=str2num(field); % slice number
      if(isempty(TP.SlicesUsed)) which=[];
      else which=find(TP.SlicesUsed==snum); end;
	if(isempty(which)) % then this is the anterior coordinate
	  TP.SlicesUsed=[TP.SlicesUsed snum];
	  TP.Slice_ant=[TP.Slice_ant;fscanf(in,'%f',[1,2])];
	else % slice already exists; then, this is the posterior coord.
	  TP.Slice_post(which,:)=fscanf(in,'%f',[1,2]);
	end
    else % not a slice; in this case it is a field
      switch(lower(field))
       case 'ac', TP.AC=fscanf(in,'%i',[1,2]);
       case 'pc', TP.PC=fscanf(in,'%i',[1,2]);
       case 'top_of_vca', TP.VAC_sup=fscanf(in,'%i',[1,2]);
       case 'top_of_vpa', TP.VPC_sup=fscanf(in,'%i',[1,2]);
       case 'ant_of_acpc', TP.ACPC_ant=fscanf(in,'%i',[1,2]);
       case 'post_of_acpc', TP.ACPC_post=fscanf(in,'%i',[1,2]);
       case 'hac', TP.HAC=fscanf(in,'%i',[1,2]);
       case 'hpc', TP.HPC=fscanf(in,'%i',[1,2]);
       case 'left_of_hca', TP.HAC_left=fscanf(in,'%i',[1,2]);
       case 'right_of_hca', TP.HAC_right=fscanf(in,'%i',[1,2]);
       case 'left_of_hcp', TP.HPC_left=fscanf(in,'%i',[1,2]);
       case 'right_of_hcp', TP.HPC_right=fscanf(in,'%i',[1,2]);
       case 'inf_trace_ac', TP.inf_trace_AC=fscanf(in,'%i',[1,2]);
       case 'inf_trace_pc', TP.inf_trace_PC=fscanf(in,'%i',[1,2]);
       case 'inf_trace_scale', TP.inf_trace_scale=fscanf(in,'%f',1);
       case 'l_trace_ac', TP.L_sag_trace_AC=fscanf(in,'%i',[1,2]);
       case 'l_trace_pc', TP.L_sag_trace_PC=fscanf(in,'%i',[1,2]);
       case 'r_trace_ac', TP.R_sag_trace_AC=fscanf(in,'%i',[1,2]);
       case 'r_trace_pc', TP.R_sag_trace_PC=fscanf(in,'%i',[1,2]);
       case 'sag_trace_scale', TP.sag_trace_scale=fscanf(in,'%f',1);
       otherwise
	fprintf(1,'\nError: I don''t recognise the field %s\n',field);
      end % switch
    end % not a slice
    field=fscanf(in,'%s',1);
  end % while reading through the file
  fclose(in);
  
  % Now, do a heck of a lot of error-trapping
  errtemplate='\nError: %s never assigned\n';
  
  fld='AC'; val=fld; if(~isfield(TP,fld)), fprintf(1,errtemplate,val); end;
  fld='PC'; val=fld; if(~isfield(TP,fld)), fprintf(1,errtemplate,val); end;
  fld='VAC_sup'; val='top_of_VCA never assigned';
  if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='VPC_sup'; val='top_of_VPA';
  if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='ACPC_ant'; val='ant_of_ACPC';
  if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='ACPC_post'; val='post_of_ACPC';
  if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='HAC'; val=fld; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='HPC'; val=fld; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='HAC_left'; val='left_of_HAC';
  if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='HAC_right'; val='right_of_HAC';
  if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='inf_trace_AC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='inf_trace_PC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='inf_trace_scale'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='L_sag_trace_AC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='L_sag_trace_PC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='R_sag_trace_AC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='R_sag_trace_PC'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  fld='sag_trace_scale'; if(~isfield(TP,fld)), fprintf(1,errtemplate,fld); end;
  
  if(size(TP.Slice_post,1)<length(TP.SlicesUsed))
    fprintf(1,'\nError: Slice %i has no post. coord.\n',SlicesUsed(end));
  end
  for snum=1:(length(TP.SlicesUsed-1))
    if(TP.Slice_post(snum,:)==[0 0])
      fprintf(1,'\nError: Slice %i has no post. coord.\n',TP.SlicesUsed(snum));
    end
  end
end % Doing the slice-wise method
