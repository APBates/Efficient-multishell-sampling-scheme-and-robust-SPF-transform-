function [even_SH_coeffs] = remove_odd_degree_SH(SH_coeff, L)
%this function returns the SH coeffcients with just even orders
% L is the maximum SH degree
%author: Alice Bates
%date last changed: 9/11/17

even_SH_coeffs = zeros(1,(L+1)*(L+2)*0.5);

index =1;
for m=0:L
    if mod(m,2)==0
       even_SH_coeffs(index:index+2*m)=SH_coeff(m^2+1:(m+1)^2);
        index = index + 2*m+1;
    end
end