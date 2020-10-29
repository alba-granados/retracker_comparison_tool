function [total_time] = compute_total_time_FBR(folder_data)

files = dir(strcat(folder_data,'*.DBL'));
total_time = 0;
disp(strcat('Total #files ',num2str(length(files))));
for i_file=1:length(files)
    file_name=char(files(i_file).name);
    
    init_time = datenum(file_name(20:20+14),'yyyymmddTHHMMSS');
    end_time = datenum(file_name(36:36+14),'yyyymmddTHHMMSS');
    total_time = total_time+((end_time-init_time)*24*60*60);
end


end