file_nc_zeros_original='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\results\Phase_2\test\2013\changes_L1B\L1B_ISR\agulhas\S3\data\CR2_SR_1_SRA____20130102T134943_20130102T135354_20170210T142844_isd.nc';
file_nc_zeros_no_noisy_beams='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\results\2013_L1B_ISD_S3_no_noisy_beams_removal\agulhas\data\CR2_SR_1_SRA____20130102T134943_20130102T135354_20170326T124228_isd_S3_no_noisy_beams.nc';
file_nc_no_zeros='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\results\Phase_2\test\2013\changes_L1B\L1B_ISR\agulhas\S3_nozeros\data\CR2_SR_1_SRA____20130102T134943_20130102T135354_20170210T182451_isd_S3_nozeros.nc';

N_samples=128;

i2q2_meas_zeros_original=double(ncread(file_nc_zeros_original,'i2q2_meas_ku_l1b_echo_sar_ku'));
scale_factor=double(ncread(file_nc_zeros_original,'waveform_scale_factor_l1b_echo_sar_ku'));
i2q2_meas_zeros_original=i2q2_meas_zeros_original.*repmat(scale_factor.',N_samples,1);


i2q2_meas_zeros_no_noisy_beams=double(ncread(file_nc_zeros_no_noisy_beams,'i2q2_meas_ku_l1b_echo_sar_ku'));
scale_factor=double(ncread(file_nc_zeros_no_noisy_beams,'waveform_scale_factor_l1b_echo_sar_ku'));
i2q2_meas_zeros_no_noisy_beams=i2q2_meas_zeros_no_noisy_beams.*repmat(scale_factor.',N_samples,1);


i2q2_meas_zeros_no_zeros=double(ncread(file_nc_no_zeros,'i2q2_meas_ku_l1b_echo_sar_ku'));
scale_factor=double(ncread(file_nc_no_zeros,'waveform_scale_factor_l1b_echo_sar_ku'));
i2q2_meas_zeros_no_zeros=i2q2_meas_zeros_no_zeros.*repmat(scale_factor.',N_samples,1);

max_image=10*log10(max([max(i2q2_meas_zeros_no_zeros(:)),max(i2q2_meas_zeros_no_noisy_beams(:)),max(i2q2_meas_zeros_original(:))]));
min_image=max_image-40;

figure; 
imagesc(10*log10(i2q2_meas_zeros_original.')); colormap('jet'); colorbar; caxis([min_image,max_image]);
xlabel('Samples'); ylabel('Waveforms'); title('Waveforms with zeros original')

figure; 
imagesc(10*log10(i2q2_meas_zeros_no_noisy_beams.')); colormap('jet'); colorbar; caxis([min_image,max_image]);
xlabel('Samples'); ylabel('Waveforms'); title('Waveforms with zeros (no noisy beams)');

figure; 
imagesc(10*log10(i2q2_meas_zeros_no_zeros.')); colormap('jet'); colorbar; caxis([min_image,max_image]);
xlabel('Samples'); ylabel('Waveforms'); title('Waveforms no zeros');

figure; 
imagesc(10*log10(abs(i2q2_meas_zeros_original.'./i2q2_meas_zeros_no_noisy_beams.'))); colormap('jet'); colorbar; caxis([min_image,max_image]);
xlabel('Samples'); ylabel('Waveforms'); title('Waveforms no zeros (Difference original & reprocessed no noisy beams)');

figure; 
imagesc((abs(i2q2_meas_zeros_original.'-i2q2_meas_zeros_no_noisy_beams.'))); colormap('jet'); colorbar;
xlabel('Samples'); ylabel('Waveforms'); title('Waveforms no zeros (Difference original & reprocessed no noisy beams)');



