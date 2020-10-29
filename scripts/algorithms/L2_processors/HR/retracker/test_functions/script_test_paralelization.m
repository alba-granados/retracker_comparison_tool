num_pools=3;
version_matlab=version;

total_files=8;
N_samples=[10,20,30,10,20,30,50,60];

if str2double(version_matlab(end-5:end-2))>2013
    parpool(num_pools);
else
    matlabpool('open',num_pools);
end

res={};
a={[]};
b={[]};
c={};
d={};
e={};
 parfor i_file=1:total_files
     disp(i_file);
     [res{i_file},a{i_file},c{i_file},d{i_file},e{i_file},b{i_file}]=generate_rands(N_samples(i_file));
 end
 
 %close pools
 if str2double(version_matlab(end-5:end-2))>2013
     poolobj = gcp('nocreate');
     delete(poolobj);
 else
     matlabpool('close');
 end
 
 disp('End')

