%creates and stores the Pm matrices for different m
L=11;
Pm_struct = cell(1,L);

[THETA, ~] = nsht_sampling_points(L); 
for m=L-1:-2:0 
    [P_1, Sc_1] = nsht_legmat_mex(THETA, L, m);
    P_mat_1 = P_1.*10.^Sc_1;
    Pm_struct{m+1} = P_mat_1;

    %for m-1
    if m~=0
        [P_2, Sc_2] = nsht_legmat_mex(THETA, L, m-1);
        P_mat_2 = P_2.*10.^Sc_2;
        Pm_struct{m} = P_mat_2;   
    end    
end


save('Pm_shell_4.mat','Pm_struct');