function [stack,beam_angles,beam_angles_tangent,lon,lat,time,index_stack,Stack_Beam_Index,Stack_Doppler_Beam_Angle,NLook_20Hz]=...
        return_stack_GPOD(time_in,L2_file,L1BS_folder,varargin)
    
    p = inputParser;
    p.addParamValue('index_stack_in',0);
    p.parse(varargin{:});
    index_stack=p.Results.index_stack_in;
    
    %% ----------- Define output variables as empty -----------------------
    stack=[];
    lon=[];
    lat=[];
    time=[];

    if ~isempty(index_stack) && index_stack==0
        %% ------------- Read lon lat & time information whole track L2 -------
        GPOD_lat_surf=double(ncread(L2_file,'latitude_20Hz').');
        GPOD_lon_surf=double(ncread(L2_file,'longitude_20Hz').');
        GPOD_time_TAI=double(ncread(L2_file,'TAI_Time_20Hz').');
        
        
        %% ------------- Look for the surface/stack index ---------------------
        [~,index_stack]=min(abs(GPOD_time_TAI-time_in));
        
        if isempty(index_stack)
            disp('No stack found for the defined time');
            return;
        end                        
    end
    
    
    
    %% ------------- Look for the nc stack containing the interest --------
    id_acq=strsplit(L2_file,'RES_CS_LTA__SIR1SAR_FR_');
    id_acq=char(id_acq(end));
    id_acq=id_acq(1:31);
    
    files_stacks=dir([L1BS_folder '*' id_acq '*.nc']);
    if isempty(files_stacks)
        error('No valid netcdf file containing the stacks is available in the specified folder')
    end
    for i_file=1:length(files_stacks)
        dumm=strsplit(files_stacks(i_file).name,'REC_');
        dumm=strsplit(char(dumm(end)),'_to_');
        first_stack=str2num(char(dumm(1)));
        dumm=strsplit(char(dumm(end)),'.nc');
        last_stack=str2num(char(dumm(1)));
        if index_stack<=last_stack && index_stack>=first_stack
            L1BS_file=[L1BS_folder char(files_stacks(i_file).name)];
            break;
        end
    end
    clear dumm;
    %% ------------ Extract the info from stack ---------------------------
    stacks=ncread(L1BS_file,'SAR_Stack_Data');
    stack=squeeze(stacks(:,:,index_stack-first_stack+1));
    clear stacks
    
    beam_angles_stacks=ncread(L1BS_file,'Beam_Angle_wrt_Tangent_Direction');
    beam_angles=(beam_angles_stacks(:,index_stack-first_stack+1)).';
    clear beam_angles_stacks;
    
    beam_angles_tangent_stacks=ncread(L1BS_file,'Beam_Angle_wrt_Velocity_Direction');
    beam_angles_tangent=(beam_angles_tangent_stacks(:,index_stack-first_stack+1)).';
    clear beam_angles_tangent_stacks;
    
    dumm=ncread(L1BS_file,'longitude');
    lon=(dumm(index_stack-first_stack+1));
    clear dumm;
    
    dumm=ncread(L1BS_file,'latitude');
    lat=(dumm(index_stack-first_stack+1));
    clear dumm;
    
    dumm=ncread(L1BS_file,'TAI_Time');
    time=(dumm(index_stack-first_stack+1));
    clear dumm;

    %% ------------- Extract info L2 --------------------------------------
    dumm=ncread(L2_file,'Stack_Beam_Index');
    Stack_Beam_Index=dumm(:,index_stack).';
    clear dumm;
    
    dumm=ncread(L2_file,'Stack_Doppler_Beam_Angle');
    Stack_Doppler_Beam_Angle=dumm(:,index_stack).';
    clear dumm;
    
    dumm=ncread(L2_file,'NLook_20Hz');
    NLook_20Hz=dumm(index_stack);
    clear dumm;
    
    
%     ncid = netcdf.open(L1BS_file,'NOWRITE');
%     [~, N_samples] = netcdf.inqDim(ncid,0);
%     [~, N_beams] = netcdf.inqDim(ncid,1);
%         
%     var_SAR_Stack_Data = netcdf.inqVarID(ncid,'SAR_Stack_Data');
%     SAR_Stack_Data_scale_factor=ncreadatt(L1BS_file,'SAR_Stack_Data','scale_factor');
%     SAR_Stack_Data_offset=ncreadatt(L1BS_file,'SAR_Stack_Data','add_offset');
%     
%     stack=double(netcdf.getVar(ncid,var_SAR_Stack_Data,[0 0 index_stack-1],[N_beams N_samples 1])).*SAR_Stack_Data_scale_factor+SAR_Stack_Data_offset;
%     
    
    
end