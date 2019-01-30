%run this at initalisation
%add paths
addpath('../');
addpath('support_functions');
addpath('theta_locations'); %for theta location of samples in the proposed sampling grid
addpath(genpath('process_real_data')); %contains scripts for processing real data
addpath(genpath('syn_eval_mat_file_generating_scripts')); %contains scripts for generating results for synthetic evaluation part of paper
addpath('synthetic_results'); %contains .mat files for plotting synthetic evaluation figures for paper

%install SHT 
install_sht;


