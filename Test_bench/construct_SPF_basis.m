function [basis] = construct_SPF_basis(Lmax, N, XYZB, diff_time,zeta,REAL)
% construct_SPF_basis.m returns the SPF matrix of SPF functions evaluated at sample locations
%
% inputs:
%
% Lmax       the SH band-limit of the SPF basis, max SH degree, must be even
% N           the radial band-limit of the SPF basis, max radial degree
% XYZB        the sampling Cartesian coordinates
% diff_time   the diffusion time parameter of the sampling protocol
% zeta        the scaling factor of the SPF basis
% REAL        if REAL=1 the returns the REAL SPF basis (SH real) otherwise returns the complex SPF basis (complex SH)
%
% output:
% basis       rows are sampled SPF functions. Size number of SPF functions x number of
%             samples.
%
% author: Alice Bates
% date last altered: 13/02/18
%%

%find the polar coordinates of samples
[x_phi,x_theta,~] = cart2sph(XYZB(:,1),XYZB(:,2),XYZB(:,3));

x_theta = pi/2-x_theta;
x_rad = transpose((XYZB(:,4)*1e6/(4*pi^2*diff_time*zeta)));

Rn = zeros((N+1)*(Lmax+1)*(Lmax+2)/2,length(x_theta));
%expand in radial direction first
for n=0:N
    constant_n = ( (2* factorial(n))/(zeta^1.5*gamma(n+1.5)) )^0.5;
    row = constant_n.*exp(-x_rad./2).*laguerreL(n,0.5,x_rad); 
    Rn(n*(Lmax+1)*(Lmax+2)*0.5+1:(n+1)*(Lmax+1)*(Lmax+2)*0.5,:) = repmat(row,(Lmax+1)*(Lmax+2)*0.5,1);     
end
 %then expand in the angular direction
if REAL
    [Ymat] = construct_SH_basis(Lmax, XYZB(:,1:3), 2, 'real');
    Ymat = transpose(Ymat);
else
    
    Ymat = zeros((Lmax+1)*(Lmax+2)/2,length(x_theta));
    %only even degree SH
    for l=0:2:Lmax
       for m=-l:l
           i=(l+1)*(l+2)/2-l+m;
           Ymat(i,:)=sphHarm(l,m,x_theta,x_phi);
       end
    end
end
Ymat = repmat(Ymat,(N+1),1);
basis = Ymat.*Rn;