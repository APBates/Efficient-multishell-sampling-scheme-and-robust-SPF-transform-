%Need to redo all code with errorbar - check works for SNR=3 to begin with

%plotting_synthetic_evaluation.m plots the mean erros in the spectral (NRMSE) and spatial domains, as well as the std in these errors verse crossing angle 

%% SNR = 10
x = 0.5:0.1:1;

%load .mat files

%proposed regularisation
load('SNR10_syn_eval_1_fibre_reg_proposed.mat');
reg_spectral_NRMSE_vec = spectral_NRMSE_vec;
reg_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec;
reg_spatial_NRMSE_vec = spatial_NRMSE_vec;
reg_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec;
clear spectral_NRMSE_vec spatial_NRMSE_vec std_spectral_NRMSE_vec std_spatial_NRMSE_vec

%LS regularisation
load('SNR10_syn_eval_1_fibre_reg_LS.mat');
LS_reg_spectral_NRMSE_vec = LS_spectral_NRMSE_vec;
LS_reg_std_spectral_NRMSE_vec = LS_std_spectral_NRMSE_vec;
LS_reg_spatial_NRMSE_vec = LS_spatial_NRMSE_vec;
LS_reg_std_spatial_NRMSE_vec = LS_std_spatial_NRMSE_vec;
clear LS_spectral_NRMSE_vec LS_spatial_NRMSE_vec LS_std_spectral_NRMSE_vec LS_std_spatial_NRMSE_vec

%proposed denoised
load('SNR10_syn_eval_1_fibre_denoised_proposed.mat');
denoised_spectral_NRMSE_vec = spectral_NRMSE_vec;
denoised_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec;
denoised_spatial_NRMSE_vec = spatial_NRMSE_vec;
denoised_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec;
clear spectral_NRMSE_vec spatial_NRMSE_vec std_spectral_NRMSE_vec std_spatial_NRMSE_vec

%LS denoised
load('SNR10_syn_eval_1_fibre_denoised_LS.mat');
LS_denoised_spectral_NRMSE_vec = spectral_NRMSE_vec_LS;
LS_denoised_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec_LS;
LS_denoised_spatial_NRMSE_vec = spatial_NRMSE_vec_LS;
LS_denoised_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec_LS;
clear  spectral_NRMSE_vec_LS std_spectral_NRMSE_vec_LS spatial_NRMSE_vec_LS std_spatial_NRMSE_vec_LS

%spectral figure
figure;
errorbar(x,LS_reg_spectral_NRMSE_vec,LS_reg_std_spectral_NRMSE_vec,'k-o','MarkerSize',8);
hold on;
grid on; box on;
errorbar(x,LS_denoised_spectral_NRMSE_vec,LS_denoised_std_spectral_NRMSE_vec,'g-s','MarkerSize',6);
errorbar(x,reg_spectral_NRMSE_vec,reg_std_spectral_NRMSE_vec,'r-*','MarkerSize',8);
errorbar(x,denoised_spectral_NRMSE_vec,denoised_std_spectral_NRMSE_vec,'b-d','MarkerSize',8);
xlabel('Fractional Anisotropy (FA)','FontSize',12);
xticks(x);
ylabel('NRMSE_e','FontSize',12);
%legend({'LS-Regularized','LS-Regularized-Denoised', 'nSPFt-Regularized','nSPFt-Regularized-Denoised'},'FontSize',12);

%spatial figure
figure;
errorbar(x,LS_reg_spatial_NRMSE_vec,LS_reg_std_spatial_NRMSE_vec,'k-o','MarkerSize',8);
hold on;
grid on; box on;
errorbar(x,LS_denoised_spatial_NRMSE_vec,LS_denoised_std_spatial_NRMSE_vec,'g-s','MarkerSize',6);
errorbar(x,reg_spatial_NRMSE_vec,reg_std_spatial_NRMSE_vec,'r-*','MarkerSize',8);
errorbar(x,denoised_spatial_NRMSE_vec,denoised_std_spatial_NRMSE_vec,'b-d','MarkerSize',8);
xlabel('Fractional Anisotropy (FA)','FontSize',12);
xticks(x);
ylabel('NRMSE_d','FontSize',12);
%legend({'LS-Regularized','LS-Regularized-Denoised', 'nSPFt-Regularized','nSPFt-Regularized-Denoised'},'FontSize',12);




%% SNR = 20
x = 0.5:0.1:1;

%load .mat files

%proposed regularisation
load('SNR20_syn_eval_1_fibre_reg_proposed.mat');
reg_spectral_NRMSE_vec = spectral_NRMSE_vec;
reg_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec;
reg_spatial_NRMSE_vec = spatial_NRMSE_vec;
reg_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec;
clear spectral_NRMSE_vec spatial_NRMSE_vec std_spectral_NRMSE_vec std_spatial_NRMSE_vec

%LS regularisation
load('SNR20_syn_eval_1_fibre_reg_LS.mat');
LS_reg_spectral_NRMSE_vec = LS_spectral_NRMSE_vec;
LS_reg_std_spectral_NRMSE_vec = LS_std_spectral_NRMSE_vec;
LS_reg_spatial_NRMSE_vec = LS_spatial_NRMSE_vec;
LS_reg_std_spatial_NRMSE_vec = LS_std_spatial_NRMSE_vec;
clear LS_spectral_NRMSE_vec LS_spatial_NRMSE_vec LS_std_spectral_NRMSE_vec LS_std_spatial_NRMSE_vec

%proposed denoised
load('SNR20_syn_eval_1_fibre_denoised_proposed.mat');
denoised_spectral_NRMSE_vec = spectral_NRMSE_vec;
denoised_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec;
denoised_spatial_NRMSE_vec = spatial_NRMSE_vec;
denoised_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec;
clear spectral_NRMSE_vec spatial_NRMSE_vec std_spectral_NRMSE_vec std_spatial_NRMSE_vec

%LS denoised
load('SNR20_syn_eval_1_fibre_denoised_LS.mat');
LS_denoised_spectral_NRMSE_vec = spectral_NRMSE_vec_LS;
LS_denoised_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec_LS;
LS_denoised_spatial_NRMSE_vec = spatial_NRMSE_vec_LS;
LS_denoised_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec_LS;
clear  spectral_NRMSE_vec_LS std_spectral_NRMSE_vec_LS spatial_NRMSE_vec_LS std_spatial_NRMSE_vec_LS

%spectral figure
figure;
errorbar(x,LS_reg_spectral_NRMSE_vec,LS_reg_std_spectral_NRMSE_vec,'k-o','MarkerSize',8);
hold on;
grid on; box on;
errorbar(x,LS_denoised_spectral_NRMSE_vec,LS_denoised_std_spectral_NRMSE_vec,'g-s','MarkerSize',6);
errorbar(x,reg_spectral_NRMSE_vec,reg_std_spectral_NRMSE_vec,'r-*','MarkerSize',8);
errorbar(x,denoised_spectral_NRMSE_vec,denoised_std_spectral_NRMSE_vec,'b-d','MarkerSize',8);
xlabel('Fractional Anisotropy (FA)','FontSize',12);
xticks(x);
ylabel('NRMSE_e','FontSize',12);
%legend({'LS-Regularized','LS-Regularized-Denoised', 'nSPFt-Regularized','nSPFt-Regularized-Denoised'},'FontSize',12);

%spatial figure
figure;
errorbar(x,LS_reg_spatial_NRMSE_vec,LS_reg_std_spatial_NRMSE_vec,'k-o','MarkerSize',8);
hold on;
grid on; box on;
errorbar(x,LS_denoised_spatial_NRMSE_vec,LS_denoised_std_spatial_NRMSE_vec,'g-s','MarkerSize',6);
errorbar(x,reg_spatial_NRMSE_vec,reg_std_spatial_NRMSE_vec,'r-*','MarkerSize',8);
errorbar(x,denoised_spatial_NRMSE_vec,denoised_std_spatial_NRMSE_vec,'b-d','MarkerSize',8);
xlabel('Fractional Anisotropy (FA)','FontSize',12);
xticks(x);
ylabel('NRMSE_d','FontSize',12);
%legend({'LS-Regularized','LS-Regularized-Denoised', 'nSPFt-Regularized','nSPFt-Regularized-Denoised'},'FontSize',12);



%% SNR = 30
x = 0.5:0.1:1;

%load .mat files

%proposed regularisation
load('SNR30_syn_eval_1_fibre_reg_proposed.mat');
reg_spectral_NRMSE_vec = spectral_NRMSE_vec;
reg_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec;
reg_spatial_NRMSE_vec = spatial_NRMSE_vec;
reg_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec;
clear spectral_NRMSE_vec spatial_NRMSE_vec std_spectral_NRMSE_vec std_spatial_NRMSE_vec

%LS regularisation
load('SNR30_syn_eval_1_fibre_reg_LS.mat');
LS_reg_spectral_NRMSE_vec = LS_spectral_NRMSE_vec;
LS_reg_std_spectral_NRMSE_vec = LS_std_spectral_NRMSE_vec;
LS_reg_spatial_NRMSE_vec = LS_spatial_NRMSE_vec;
LS_reg_std_spatial_NRMSE_vec = LS_std_spatial_NRMSE_vec;
clear LS_spectral_NRMSE_vec LS_spatial_NRMSE_vec LS_std_spectral_NRMSE_vec LS_std_spatial_NRMSE_vec

%proposed denoised
load('SNR30_syn_eval_1_fibre_denoised_proposed.mat');
denoised_spectral_NRMSE_vec = spectral_NRMSE_vec;
denoised_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec;
denoised_spatial_NRMSE_vec = spatial_NRMSE_vec;
denoised_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec;
clear spectral_NRMSE_vec spatial_NRMSE_vec std_spectral_NRMSE_vec std_spatial_NRMSE_vec

%LS denoised
load('SNR30_syn_eval_1_fibre_denoised_LS.mat');
LS_denoised_spectral_NRMSE_vec = spectral_NRMSE_vec_LS;
LS_denoised_std_spectral_NRMSE_vec = std_spectral_NRMSE_vec_LS;
LS_denoised_spatial_NRMSE_vec = spatial_NRMSE_vec_LS;
LS_denoised_std_spatial_NRMSE_vec = std_spatial_NRMSE_vec_LS;
clear  spectral_NRMSE_vec_LS std_spectral_NRMSE_vec_LS spatial_NRMSE_vec_LS std_spatial_NRMSE_vec_LS

%spectral figure
figure;
errorbar(x,LS_reg_spectral_NRMSE_vec,LS_reg_std_spectral_NRMSE_vec,'k-o','MarkerSize',8);
hold on;
grid on; box on;
errorbar(x,LS_denoised_spectral_NRMSE_vec,LS_denoised_std_spectral_NRMSE_vec,'g-s','MarkerSize',6);
errorbar(x,reg_spectral_NRMSE_vec,reg_std_spectral_NRMSE_vec,'r-*','MarkerSize',8);
errorbar(x,denoised_spectral_NRMSE_vec,denoised_std_spectral_NRMSE_vec,'b-d','MarkerSize',8);
xlabel('Fractional Anisotropy (FA)','FontSize',12);
xticks(x);
ylabel('NRMSE_e','FontSize',12);
%legend({'LS-Regularized','LS-Regularized-Denoised', 'nSPFt-Regularized','nSPFt-Regularized-Denoised'},'FontSize',12);

%spatial figure
figure;
errorbar(x,LS_reg_spatial_NRMSE_vec,LS_reg_std_spatial_NRMSE_vec,'k-o','MarkerSize',8);
hold on;
grid on; box on;
errorbar(x,LS_denoised_spatial_NRMSE_vec,LS_denoised_std_spatial_NRMSE_vec,'g-s','MarkerSize',6);
errorbar(x,reg_spatial_NRMSE_vec,reg_std_spatial_NRMSE_vec,'r-*','MarkerSize',8);
errorbar(x,denoised_spatial_NRMSE_vec,denoised_std_spatial_NRMSE_vec,'b-d','MarkerSize',8);
xlabel('Fractional Anisotropy (FA)','FontSize',12);
xticks(x);
ylabel('NRMSE_d','FontSize',12);
%legend({'LS-Regularized','LS-Regularized-Denoised', 'nSPFt-Regularized','nSPFt-Regularized-Denoised'},'FontSize',12);

