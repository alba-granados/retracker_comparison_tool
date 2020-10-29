version_matlab=version;
num_pools=1;

ftpPath     = 'science-pds.cryosat.esa.int';
user        = 'cryosat081';
pswrd       = '4MsIuLNr';


type_product_out = 'SIR_SAR_L2';
type_product_in = 'SIR_SAR_L2';
type_input_file='nc';
mode='SAR';


inputDir_vec    = {'C:/Users/eduard.makhoul/Desktop/temporal/SCOOP/'};
outputDir_vec   = {'C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L2_ESA/north_sea/2013/'};


n_folders=length(inputDir_vec);


if num_pools~=1
    %create pools
    if str2double(version_matlab(end-5:end-2))>2013
        parpool(num_pools);    
    else
        matlabpool('open',num_pools);
    end
    %% ------------- Loop per folder ------------------------
    parfor i_folder=1:n_folders
        try
            download_CR2_ftp (ftpPath,user,pswrd,type_product_out,type_product_in,char(inputDir_vec(i_folder)),char(outputDir_vec(i_folder)),type_input_file,mode)
        catch
            continue;
        end
    end
    %close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end    
else
    for i_folder=1:n_folders
        try
            download_CR2_ftp (ftpPath,user,pswrd,type_product_out,type_product_in,char(inputDir_vec(i_folder)),char(outputDir_vec(i_folder)),type_input_file,mode)
        catch
            continue;
        end
    end
end

