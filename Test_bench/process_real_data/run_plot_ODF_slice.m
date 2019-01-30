% Script to plot the fODF in a brain slice
clear; clc
path_to_mask = 'process_real_data/opt_dim_multi_mask.nii';

folder  = {'process_real_data/ODF_folder/'};

%_________________________________________________________________________%
for i=1:length(folder)
    folder_i = char(folder(i));
    display (['Creating figure for path ' num2str(i)]);
    display (['Folder -> ' folder_i]);
    display (' ');
    % ---------------------------------------------------------------------
    path_to_ODF =    [folder_i '/ODF_matrix_4D.nii'];    
    path_to_GFA =    [folder_i '/dt_opt_dim_multi_fa.nii'];

    % ---------------------------------------------------------------------
    %type = 'y'; slice = 55;
    type = 'y'; slice = 72;%73;
    %type = 'y'; slice = 56;
% 
     x1 = 45; x2 = 70;
     y1 = 53; y2 = 70;

%     x1 = 65; x2 = 82;
%     y1 = 54; y2 = 60;


    plot_ODFs_SF(path_to_ODF, path_to_mask, path_to_GFA, slice, type);

    axis([x1 x2 slice-2 slice+2 y1 y2]); box on; axis off
    line([x1 x2],[slice slice],[y1 y1],'LineWidth',12,'Color',[.2 .2 .8]')
    line([x1 x2],[slice slice],[y2 y2],'LineWidth',12,'Color',[.2 .2 .8]')
    line([x1 x1],[slice slice],[y1 y2],'LineWidth',12,'Color',[.2 .2 .8]')
    line([x2 x2],[slice slice],[y1 y2],'LineWidth',12,'Color',[.2 .2 .8]')
    figure_name = [folder_i '/plot_ODF_roi_slice_' type '_' type num2str(slice)];
    
    % If you install "export_fig" from matlab exchange, then you can use the
    % following line to export your figure:
    % export_fig (figure_name, '-jpg', '-nocrop', '-transparent', '-a1', '-r300', '-nocrop');
    % close all;
    % ---------------------------------------------------------------------
end