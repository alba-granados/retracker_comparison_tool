function []=change_file_name_numerical(input_path,output_path,file_ext,token_left,token_right,numeric_string_fromat)


filesBulk.inputPath       =   input_path;
mkdir(output_path);

filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
filesBulk.indexFilesNC=find(~cellfun(@isempty,strfind(aux,file_ext)));
filesBulk.nFilesNC=length(filesBulk.indexFilesNC);
filesBulk.NCFiles=filesBulk.inputFiles(filesBulk.indexFilesNC);

for i_file=1:filesBulk.nFilesNC
    filename=char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name);
    
    filtered_name = strsplit(filename,token_left);
    if length(filtered_name)==1
        continue;
    end
    filtered_name = strsplit(char(filtered_name(end)),token_right);
    if length(filtered_name)==1
        continue;
    end
    
    
    numerical_value = str2num(char(filtered_name(1)));
    
    filename_renamed = strrep(filename,strcat(token_left,num2str(numerical_value)),...
                            strcat(token_left,num2str(numerical_value,numeric_string_fromat)));
    
    
    
    %copy this file in the output folder
    movefile(strcat(input_path,filename),strcat(output_path,filename_renamed));
    
end



end