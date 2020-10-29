folder1    = '\\X8DTL\data\DATA\SCOOP\copy_disk_ESA_2013\agulhas\';
folder2    = '\\X8DTL\data\DATA\SCOOP\GPOD\l2\agulhas\2013_SAMOSA2_ADAPTIVE_PTR\';

file_ext_folder1 = '.DBL';
file_ext_folder2 = '.nc';

inputFiles=dir([folder1 '*' file_ext_folder1]);


for iFile=1:length(inputFiles)
    string_data = inputFiles(iFile).name(20:50);
    inputFile2 = dir([folder2 '*' string_data '*' file_ext_folder2]);
    if isempty(inputFile2)
        disp(inputFiles(iFile).name)
    end
    
end