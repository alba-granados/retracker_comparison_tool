function []=move_file_token(input_path,output_path,file_ext,token)


filesBulk.inputPath       =   input_path;
mkdir(output_path);

filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
filesBulk.indexFilesNC=find(~cellfun(@isempty,strfind(aux,[file_ext])));
filesBulk.nFilesNC=length(filesBulk.indexFilesNC);
filesBulk.NCFiles=filesBulk.inputFiles(filesBulk.indexFilesNC);

for i_file=1:filesBulk.nFilesNC
    filename=char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name);

    filtered_name = regexp(filename,token, 'once');
    
    if ~isempty(filtered_name)
        movefile(strcat(input_path,filename),strcat(output_path,filename));        
    end
%     %copy this file in the output folder

    
end



end