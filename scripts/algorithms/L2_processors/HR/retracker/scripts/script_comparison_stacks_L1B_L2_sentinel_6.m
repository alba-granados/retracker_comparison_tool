filename_stacks_L1B = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_CPP\lat0\l1b_hr_product_lat00.mat';
filename_stacks_modelled_L2 = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\results\L2\lat0\MAT\data\l1b_hr_product_lat00_Stacks_info.mat';

path_results='C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\results\L2\lat0\MAT\plots\fitted_stacks\';
mkdir(path_results);


figure_format='jpg';
res_fig='-r300';

switch lower(figure_format)
    case 'pdf'
        file_ext='.pdf';
        print_file='-dpdf';
    case 'eps'
        file_ext='.eps';
        print_file='-depsc';
    case 'png'
        file_ext='.png';
        print_file='-dpng';
    case 'jpg'
        file_ext='.jpg';
        print_file='-djpeg';        
end

stacks_L1B=load(filename_stacks_L1B);
stacks_modelled_L2=load(filename_stacks_modelled_L2);


set_default_plot;

s = size(stacks_L1B.beams_masked);
N_surfaces = s(1);
N_max_beams = s(1);
%% ------------ Stacks -----------------------------------------
max_image = 0.0;%10*log10(max(stacks_L1B.beams_masked(:)));
min_image = -20.0; %max_image-40;
for i_surf=1:N_surfaces
    max_stack_L1B = max(max(stacks_L1B.beams_masked(i_surf,:,:)));
    max_stack_L2 = max(max(stacks_modelled_L2.modelled_stacks(i_surf,:,:)));
    f1=figure;
    subplot(1,2,1)
    surfc(10*log10(squeeze(stacks_L1B.beams_masked(i_surf,:,:)./max_stack_L1B)),...
        'FaceColor','texturemap','EdgeColor','none'); view(2); grid off;
    colormap('jet'); c=colorbar; ylabel(c,'[dB]');
    caxis([min_image,max_image]);
    title('L1B-S','Interpreter','none')
    xlabel('Samples'); ylabel('Beams');
    
    subplot(1,2,2)
    surfc(10*log10(squeeze(stacks_modelled_L2.modelled_stacks(i_surf,:,:)./max_stack_L2)),...
        'FaceColor','texturemap','EdgeColor','none'); view(2); grid off;
    colormap('jet'); c=colorbar; ylabel(c,'[dB]');
    caxis([min_image,max_image]);
    title('Modelled L2','Interpreter','none')
    xlabel('Surfaces'); ylabel('Beams');
    [axT,hT]=suplabel('Stacks comparison','t');
    print(print_file,res_fig,strcat(path_results,'stack_comparison_surf_',...
        num2str(i_surf,'%04.0f'),file_ext));
    close(f1)
end

