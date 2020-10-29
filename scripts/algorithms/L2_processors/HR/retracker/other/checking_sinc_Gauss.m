function checking_sinc_Gauss
% set fitting option
close all;
clear all;
file_ext='.eps';
print_file='-depsc';
res_fig='-r300';
path_results='C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\documentation\WPs\WP5000\ATBD_retracker\Review_retracker\figures\'
normalized_window_active =1;
options     =   optimset('Algorithm', 'levenberg-marquardt','Display','off');
fit_params_init_onesided=[0.75,0.25]; %amplitude and
%fit_params_init=[1.0,64*64/2,0.5]; %amplitude and
%% -------------- PTR description as per Chris paper ----------------------
% A one-sided Gaussian/PTR description
N=64; %number of pulses 
%N_bis=512;
m=(1-N/2):1:N/2;
%m_bis=(1-N_bis/2):1:N_bis/2;
%w_boxcar_bis=ones(1,N_bis);
%Gamma^2(xsi)=A_gexp(-xsi^2/(2*sigma_g^2))
xsi=0:1e-2:10; %range as defined by Chris in Fig.4 0 to 2; frequency components (oversampling)
w_boxcar=ones(1,N);
%w_hamming=0.54-0.46*cos(2*pi*(m-N/2)/(N));
%w_hanning=0.5*(1.0-cos(2*pi*(m-N/2)/(N)));
w_hamming=hamming(N,'periodic').';
w_hanning=hann(N,'periodic').';
w_hamming_sym=hamming(N,'symmetric').';
w_hanning_sym=hann(N,'symmetric').';


gamma_boxcar=zeros(1,length(xsi));
gamma_hamming=zeros(1,length(xsi));
gamma_hanning=zeros(1,length(xsi));
gamma_hamming_sym=zeros(1,length(xsi));
gamma_hanning_sym=zeros(1,length(xsi));
for i_xsi=1:length(xsi)
    gamma_boxcar(i_xsi)=abs(1/N*sum(w_boxcar.*exp(1i*2*pi*xsi(i_xsi)*m/N)));
    %gamma_boxcar_bis(i_xsi)=abs(1/N_bis*sum(w_boxcar_bis.*exp(1i*2*pi*xsi(i_xsi)*m_bis/N_bis)));
    gamma_hamming(i_xsi)=abs(1/N*sum(w_hamming.*exp(1i*2*pi*xsi(i_xsi)*m/N)));    
    gamma_hanning(i_xsi)=abs(1/N*sum(w_hanning.*exp(1i*2*pi*xsi(i_xsi)*m/N)));    
    gamma_hamming_sym(i_xsi)=abs(1/N*sum(w_hamming_sym.*exp(1i*2*pi*xsi(i_xsi)*m/N)));    
    gamma_hanning_sym(i_xsi)=abs(1/N*sum(w_hanning_sym.*exp(1i*2*pi*xsi(i_xsi)*m/N))); 
end

if normalized_window_active
    gamma_hamming=gamma_hamming./(mean(w_hamming));
    gamma_hanning=gamma_hanning./(mean(w_hanning));
    gamma_hamming_sym=gamma_hamming_sym./(mean(w_hamming_sym));
    gamma_hanning_sym=gamma_hanning_sym./(mean(w_hanning_sym));
end

% figure;
% title('One-sided PTR comparison')
% plot(xsi,gamma_boxcar.^2,'-b');
% hold on;
% plot(xsi,gamma_hamming.^2,'--r');
% plot(xsi,gamma_hanning.^2,'.-g');
% legend('Boxcar','Hamming','Hanning')
% 
% figure;
% title('One-sided PTR comparison')
% plot(xsi,20*log10(gamma_boxcar),'-b');
% axis([0,10,-60,0]);
% hold on;
% plot(xsi,20*log10(gamma_hamming),'--r');
% plot(xsi,20*log10(gamma_hamming_sym),'.-g');
% legend('Boxcar','Periodic','Symmetric')
% ylabel('[dB]')

% figure;
% title('One-sided PTR comparison')
% plot(xsi,20*log10(gamma_boxcar),'-b');
% axis([0,10,-60,0]);
% hold on;
% plot(xsi,20*log10(gamma_hanning),'--r');
% plot(xsi,20*log10(gamma_hanning_sym),'.-g');
% legend('Boxcar','Periodic','Symmetric')
% ylabel('[dB]')


figure;
plot(xsi,(gamma_boxcar.^2),'-b');
ylim([0,1]);
title('One-sided PTR comparison')
xlabel('\xi');
ylabel('Power')
hold on;
plot(xsi,(gamma_hamming.^2),'-r');
plot(xsi(1:10:end),(gamma_hamming_sym(1:10:end).^2),'or');
plot(xsi,(gamma_hanning).^2,'-g');
plot(xsi(1:10:end),(gamma_hanning_sym(1:10:end).^2),'*g');
legend('No window','Hamming Periodic','Hamming Symmetric','Hanning Periodic','Hanning Symmetric')




figure;
plot(xsi,20*log10(gamma_boxcar),'-b');
ylim([-60,0]);
title('One-sided PTR comparison')
xlabel('\xi');
ylabel('[dB]')
hold on;
plot(xsi,20*log10(gamma_hamming),'-r');
plot(xsi(1:10:end),20*log10(gamma_hamming_sym(1:10:end)),'or');
plot(xsi,20*log10(gamma_hanning),'-g');
plot(xsi(1:10:end),20*log10(gamma_hanning_sym(1:10:end)),'*g');
legend('No window','Hamming Periodic','Hamming Symmetric','Hanning Periodic','Hanning Symmetric')
print(print_file,res_fig,strcat(path_results,'one_sided_PTRs',file_ext))

%fitting
function Gauss_func=gaussian_one_sided(params_g,x)
    Gauss_func=params_g(1)*exp(-(x).^2/(2.0*(params_g(2)).^2));
end

%fitting to a Gaussian
%lsq
mpfun=@(params_g,x)gaussian_one_sided(params_g,x);

%% -------------------------- No window --------------------------------------
%--------------------------------------------------------------------------
%lsq minimization
params_g_boxcar_lsq =   lsqcurvefit (mpfun,fit_params_init_onesided,xsi,...
                    gamma_boxcar.^2,[ ],[ ],options);                                                
Gauss_approx_boxcar_lsq=gaussian_one_sided(params_g_boxcar_lsq,xsi);                
rmse_gauss_boxcar_lsq=sqrt(mean((gamma_boxcar.^2-Gauss_approx_boxcar_lsq).^2));

%fmin minimization
fminfun=@(fit_params)sum((gaussian_one_sided(fit_params,xsi)-gamma_boxcar.^2).^2);
params_g_boxcar_fmin=fminsearchbnd (fminfun,fit_params_init_onesided);
Gauss_approx_boxcar_fmin=gaussian_one_sided(params_g_boxcar_fmin,xsi);                
rmse_gauss_boxcar_fmin=sqrt(mean((gamma_boxcar.^2-Gauss_approx_boxcar_fmin).^2));

% %minimization with constrains (same 3dB width as boxcar)
% function Gauss_func=gaussian_one_sided_const(params_g,x,x_ref,y_ref)
%     Gauss_func=params_g(1)*exp(-(x).^2/(2.0*(params_g(2)).^2));
%     for i_ref=1:length(x_ref)
%         ref_idx=find(x==x_ref(i_ref), 1);
%         if ~isempty(ref_idx)
%             Gauss_func(ref_idx)=y_ref(i_ref);
%         end
%     end
% end
% idx_3dB=find(10*log10(gamma_boxcar.^2) <= (max(10*log10(gamma_boxcar.^2))-3.0));
% b_const=gamma_boxcar(idx_3dB(1)).^2;
% fminfun=@(fit_params)sum((gaussian_one_sided_const(fit_params,xsi,[0,idx_3dB(1)],[1,b_const])-gamma_boxcar.^2).^2);
% params_g_boxcar_fmin_const=fmincon(fminfun,fit_params_init_onesided);
% Gauss_approx_boxcar_fmin_const=gaussian_one_sided(params_g_boxcar_fmin_const,xsi);                
% rmse_gauss_boxcar_fmin_const=sqrt(mean((gamma_boxcar.^2-Gauss_approx_boxcar_fmin_const).^2));

figure; 
plot(xsi,(gamma_boxcar.^2),'-b');
hold on; plot(xsi,Gauss_approx_boxcar_lsq,'--r');
% plot(xsi,Gauss_approx_boxcar_fmin,'.-g');
% plot(xsi,Gauss_approx_boxcar_fmin_const,':m');
title('PTR power gaussian approximation (No window)')
xlabel('\xi'); ylabel('Power');
text(0.5,0.85*max(gamma_boxcar.^2),strcat('RMSE (lsq) = ',num2str(rmse_gauss_boxcar_lsq)))
text(1.5,0.55*max(gamma_boxcar.^2),{strcat('A (lsq) = ',num2str(params_g_boxcar_lsq(1))),strcat('\sigma (lsq) = ',num2str(params_g_boxcar_lsq(2)))},'Interpreter','Tex')
%text(0.5,0.75*max(gamma_boxcar.^2),strcat('RMSE (fmin) = ',num2str(rmse_gauss_boxcar_fmin)))
% text(1.5,0.35*max(gamma_boxcar.^2),{strcat('A (fmin) = ',num2str(params_g_boxcar_fmin(1))),strcat('\sigma (fmin) = ',num2str(params_g_boxcar_fmin(2)))},'Interpreter','Tex')
% text(0.5,0.65*max(gamma_boxcar.^2),strcat('RMSE (fmin-const) = ',num2str(rmse_gauss_boxcar_fmin_const)))
% text(1.5,0.15*max(gamma_boxcar.^2),{strcat('A (fmin-const) = ',num2str(params_g_boxcar_fmin_const(1))),strcat('\sigma (fmin-const) = ',num2str(params_g_boxcar_fmin_const(2)))},'Interpreter','Tex')
%legend('PTR boxcar','Gaussian approx (lsq)','Gaussian approx (fmin)','Gaussian approx (fmin-constrains 3dB)');
legend('PTR boxcar','Gaussian approx (lsq)');
print(print_file,res_fig,strcat(path_results,'PTR_gauss_fit_no_window',file_ext))

figure; 
plot(xsi,Gauss_approx_boxcar_lsq-gamma_boxcar.^2,'--r');
hold on;
ylim([-0.06,0.06]);
%plot(xsi,Gauss_approx_boxcar_fmin-gamma_boxcar.^2,'.-g');
% plot(xsi,Gauss_approx_boxcar_fmin_const-gamma_boxcar.^2,':m');
title('Error caused by Gaussian approximation (No window)')
xlabel('\xi'); ylabel('Power');
%legend('Gaussian approx (lsq)','Gaussian approx (fmin)','Gaussian approx (fmin-constrains 3dB)');
%legend('PTR boxcar','Gaussian approx (lsq)');
print(print_file,res_fig,strcat(path_results,'PTR_gauss_fit_no_window_error',file_ext))



%% ------------------------ Hamming ---------------------------------------
%--------------------------------------------------------------------------
%lsq
params_g_hamming_lsq =   lsqcurvefit (mpfun,fit_params_init_onesided,xsi,...
                    gamma_hamming.^2,[ ],[ ],options);
                
Gauss_approx_hamm_lsq=gaussian_one_sided(params_g_hamming_lsq,xsi);
rmse_gauss_hamm_lsq=sqrt(mean((gamma_hamming.^2-Gauss_approx_hamm_lsq).^2));

%fmin minimization
fminfun=@(fit_params)sum((gaussian_one_sided(fit_params,xsi)-gamma_hamming.^2).^2);
params_g_hamming_fmin=fminsearchbnd (fminfun,fit_params_init_onesided);
Gauss_approx_hamm_fmin=gaussian_one_sided(params_g_hamming_fmin,xsi);                
rmse_gauss_hamm_fmin=sqrt(mean((gamma_hamming.^2-Gauss_approx_hamm_fmin).^2));


figure;
plot(xsi,(gamma_hamming.^2),'-b');
hold on; plot(xsi,Gauss_approx_hamm_lsq,'--r')
% plot(xsi,Gauss_approx_hamm_fmin,'.-g');
title('PTR power gaussian approximation (Hamming window)')
xlabel('\xi'); ylabel('Power');
text(0.5,0.65*max(gamma_hamming.^2),strcat('RMSE (lsq) = ',num2str(rmse_gauss_hamm_lsq)))
text(1.5,0.3*max(gamma_hamming.^2),{strcat('A (lsq) = ',num2str(params_g_hamming_lsq(1))),strcat('\sigma (lsq) = ',num2str(params_g_hamming_lsq(2)))},'Interpreter','Tex')
% text(0.5,0.5*max(gamma_hamming.^2),strcat('RMSE (fmin) = ',num2str(rmse_gauss_hamm_fmin)))
% text(1.5,0.15*max(gamma_hamming.^2),{strcat('A (fmin) = ',num2str(params_g_hamming_fmin(1))),strcat('\sigma (fmin) = ',num2str(params_g_hamming_fmin(2)))},'Interpreter','Tex')
%legend('PTR Hamming','Gaussian approx (lsq)','Gaussian approx (fmin)');
legend('PTR Hamming','Gaussian approx (lsq)');
print(print_file,res_fig,strcat(path_results,'PTR_gauss_fit_Hamm',file_ext))

figure; 
plot(xsi,Gauss_approx_hamm_lsq-gamma_hamming.^2,'--r');
hold on;
ylim([-0.06,0.06])
% plot(xsi,Gauss_approx_hamm_fmin-gamma_hamming.^2,'.-g');
title('Error caused by Gaussian approximation (Hamming window)')
xlabel('\xi'); ylabel('Power');
%legend('Gaussian approx (lsq)','Gaussian approx (fmin)');
%legend('PTR Hamming','Gaussian approx (lsq)');
print(print_file,res_fig,strcat(path_results,'PTR_gauss_fit_Hamm_error',file_ext))
                              
%% ------------------------ Hanning ---------------------------------------
%--------------------------------------------------------------------------
%lsq                
params_g_hanning_lsq =   lsqcurvefit (mpfun,fit_params_init_onesided,xsi,...
                    gamma_hanning.^2,[ ],[ ],options);
Gauss_approx_hann_lsq=gaussian_one_sided(params_g_hanning_lsq,xsi);
rmse_gauss_hann_lsq=sqrt(mean((gamma_hanning.^2-Gauss_approx_hann_lsq).^2));

%fmin minimization
fminfun=@(fit_params)sum((gaussian_one_sided(fit_params,xsi)-gamma_hanning.^2).^2);
params_g_hanning_fmin=fminsearchbnd (fminfun,fit_params_init_onesided);
Gauss_approx_hann_fmin=gaussian_one_sided(params_g_hanning_fmin,xsi);                
rmse_gauss_hann_fmin=sqrt(mean((gamma_hanning.^2-Gauss_approx_hann_fmin).^2));

figure; 
plot(xsi,(gamma_hanning.^2),'-b');
hold on; plot(xsi,Gauss_approx_hann_lsq,'--r')
%plot(xsi,Gauss_approx_hann_fmin,'.-g');
title('PTR power gaussian approximation (Hanning window)')
xlabel('\xi'); ylabel('Power');
text(0.5,0.65*max(gamma_hanning.^2),strcat('RMSE (lsq) = ',num2str(rmse_gauss_hann_lsq)))
text(1.5,0.3*max(gamma_hanning.^2),{strcat('A (lsq) = ',num2str(params_g_hanning_lsq(1))),strcat('\sigma (lsq) = ',num2str(params_g_hanning_lsq(2)))},'Interpreter','Tex')
% text(0.5,0.5*max(gamma_hanning.^2),strcat('RMSE (fmin) = ',num2str(rmse_gauss_hann_fmin)))
% text(1.5,0.15*max(gamma_hanning.^2),{strcat('A (fmin) = ',num2str(params_g_hanning_fmin(1))),strcat('\sigma (fmin) = ',num2str(params_g_hanning_fmin(2)))},'Interpreter','Tex')
%legend('PTR Hanning','Gaussian approx (lsq)','Gaussian approx (fmin)');
legend('PTR Hanning','Gaussian approx (lsq)');
print(print_file,res_fig,strcat(path_results,'PTR_gauss_fit_Hann',file_ext))

figure; 
plot(xsi,Gauss_approx_hann_lsq-gamma_hanning.^2,'--r');
ylim([-0.06,0.06]);
%hold on;
%plot(xsi,Gauss_approx_hann_fmin-gamma_hanning.^2,'.-g');
title('Error caused by Gaussian approximation (Hanning window)')
xlabel('\xi'); ylabel('Power');
%legend('Gaussian approx (lsq)','Gaussian approx (fmin)');
%legend('PTR Hanning','Gaussian approx (lsq)');
print(print_file,res_fig,strcat(path_results,'PTR_gauss_fit_Hann_error',file_ext))

disp('end');
% figure; 
% plot(20*log10(gamma_boxcar),'-b');
% axis([1,N,-60,0])
% hold on;
% plot(20*log10(gamma_hamming),'--r');
% plot(20*log10(gamma_hanning),'og');
% ylabel('[dB]');

% %% -------------- PTR from FFT processing point view ----------------------
% %along-track: burst information
% %definition 
% N_pulses=64;
% burst=ones(1,N_pulses);
% hamming_w=hamming(N_pulses+1); hamming_w=hamming_w(2:end).';
% hanning_w=hanning(N_pulses+1); hanning_w=hanning_w(2:end).';
% ZP_along_track=64; %oversampling factor to allow a better visualization/fitting
% 
% PTR_boxcar=fftshift(abs(1.0/sqrt(ZP_along_track*N_pulses).*fft(burst,ZP_along_track*N_pulses)));
% PTR_hamming=fftshift(abs(1.0/sqrt(ZP_along_track*N_pulses).*fft(burst.*hamming_w,ZP_along_track*N_pulses)));
% PTR_hanning=fftshift(abs(1.0/sqrt(ZP_along_track*N_pulses).*fft(burst.*hanning_w,ZP_along_track*N_pulses)));
% 
% 
% figure; 
% title('PTR power (non-normalized windows)');
% plot((PTR_boxcar.^2),'-b');
% hold on;
% plot((PTR_hamming.^2),'--r');
% plot((PTR_hanning.^2),'.-g');
% ylabel('Power');
% legend('No window','Hamming','Hanning');
% 
%                 
% function Gauss_func=gaussian_descriptor(params_g,x)
%     Gauss_func=params_g(1)*exp(-(x-params_g(2)).^2/(2.0*(params_g(3)).^2));
% end
% 
% %fitting to a Gaussian
% x=1:1:ZP_along_track*N_pulses;
% mpfun=@(params_g,x)gaussian_descriptor(params_g,x);
% params_g_boxcar =   lsqcurvefit (mpfun,fit_params_init,x,...
%                     PTR_boxcar.^2,[ ],[ ],options);
%                 
% figure; 
% title('PTR power gaussian approximation (boxcar window)')
% plot((PTR_boxcar.^2),'-b');
% axis([1600,2500,0.0,1.05])
% Gauss_approx_boxcar=gaussian_descriptor(params_g_boxcar,x);
% hold on; plot(Gauss_approx_boxcar,'--r')
% legend('PTR boxcar','Gaussian approx');
% 
% params_g_hamming =   lsqcurvefit (mpfun,fit_params_init,x,...
%                     PTR_hamming.^2,[ ],[ ],options);
% figure; 
% title('PTR power gaussian approximation (Hamming window)')
% plot((PTR_hamming.^2),'-b');
% axis([1600,2500,0.0,1.5])
% Gauss_approx_hamm=gaussian_descriptor(params_g_hamming,x);
% hold on; plot(Gauss_approx_hamm,'--r')
% legend('PTR Hamming','Gaussian approx');
%                               
%                 
% params_g_hanning =   lsqcurvefit (mpfun,fit_params_init,x,...
%                     PTR_hanning.^2,[ ],[ ],options);
% figure; 
% title('PTR power gaussian approximation (Hanning window)')
% plot((PTR_hanning.^2),'-b');
% axis([1600,2500,0.0,1.5])
% Gauss_approx_hann=gaussian_descriptor(params_g_hanning,x);
% hold on; plot(Gauss_approx_hann,'--r')
% legend('PTR Hanning','Gaussian approx');

                


end




