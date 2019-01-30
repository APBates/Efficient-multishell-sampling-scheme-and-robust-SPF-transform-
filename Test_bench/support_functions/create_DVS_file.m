%the scaling factor of this is incorrect

N=3;
NUM_SAMPLES = 94;
BMAX = 4000;
sample_i = 1;

dvs_vec = zeros(NUM_SAMPLES,3);

for i=1:N+1
    %load shell
    load(['shell_',num2str(i),'.mat']);
    shell_bvalue = XYZB(1,4);
    shell_dvs_samples = XYZB(:,1:3);
    
    %scale so that length is a fraction of the max B-value 
    b_scale = shell_bvalue/BMAX;
    shell_dvs_samples = shell_dvs_samples*b_scale;
    
    %store in dvs_vec
    dvs_vec(sample_i:sample_i+size(shell_dvs_samples,1)-1,:) = shell_dvs_samples;
    sample_i = sample_i + size(shell_dvs_samples,1);
end


%then need to copy and paste dvs_vec into a .dvs file