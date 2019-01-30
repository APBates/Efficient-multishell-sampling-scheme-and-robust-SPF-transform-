function[E] = E_expansion(Elmn,L,N,zeta,x_rad, x_theta,x_phi)
%expansion.m performs the inverse transform in the SPF basis to calculate the SPF
%coefficients. 
% Inputs:
% Elmn       the coefficients in the SPF basis stored as a N x Lmax^2
% matrix
% L          a vector containing the SH band-limit at different shells, max SH degree
% N          the radial truncation order = number of shells. Max radial
% order = N-1
% zeta       scale factor
% x_rad         the locations in q-space to determine signal (x = q^2/zeta).
% x_theta
% x_phi         
% Output:

% E          the diffusion signal at specified q space location


%author: Alice Bates


%expansion for when want to evaluate at arbitrary locations

Lmax = max(L);
Rn = zeros(N+1,1);
%expand in radial direction first
for n=0:N
    constant_n = ( (2* factorial(n))/(zeta^1.5*gamma(n+1.5)) )^0.5;
    Rn(n+1) =  constant_n.*exp(-x_rad./2).*laguerreL(n,0.5,x_rad);     
end
%then expand in the angular direction
temp = transpose(Elmn)*Rn;
Ymat = zeros(1,(Lmax+1)^2);
for i=1:(Lmax+1)^2
    l = floor(sqrt(i-1));
    m = i-1-l*(l+1);
    Ymat(i) = sphHarm(l,m,x_theta,x_phi);
end
E = Ymat*temp;