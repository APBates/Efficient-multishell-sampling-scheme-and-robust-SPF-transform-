function[E] = iSPFt(Elmn,L,N,x,zeta)
%iSPFt performs the inverse transform in the SPF basis to the diffusion
%signal at the orginal sample locations. 
% Inputs:
% Elmn       the coefficients in the SPF basis stored as a (N+1) x (Lmax+1)^2
% matrix
% L          a vector containing the SH band-limit at different shells, the
% max SH degree
% N          the radial truncation order. Max radial order. Number shells =
% N+1
% x         Gauss-Laguerre quadrature roots 
% zeta       scale factor
% Output:

% E          the diffusion signal stored as a (N+1) x (Lmax+1)*(Lmax+2)*0.5 matrix.
% i.e. each row contains samples in one shell, rows padded with zeros at
% end for shells with less samples

%author: Alice Bates
%date last altered: 15/02/18

%Using nsht inverse tranform
% Only works if same band-limit within each shell (or if want results on grid with max L in each shell)

% Lmax = max(L);
% Eqn = zeros(N+1,0.5*(Lmax+1)*(Lmax+2));
% Rqn = zeros(N+1,N+1); 
% %for each radial function
% for n2=1:N+1
%      n=n2-1;
%      constant_n = ( (2* factorial(n))/(zeta^0.5*gamma(n+1.5)) )^0.5;  
%      Rqn(:,n2) = constant_n.*exp(-x./2).*laguerreL(n,0.5,x); 
%      Eqn(n2,:) = nsht_inverse(Elmn(n2,:),Lmax+1); %provide with sig pro L
% end
%     
% E = Rqn*Eqn;


%% calculate at individual sample locations - shell by shell
Lmax = max(L);
E = zeros(N+1,(Lmax+1)*(Lmax+2)*0.5);
%for each shell
for n=1:N+1
    [theta, phis] = nsht_sampling_points(L(n)+1); % sig pro definition band-limit L
    thetas = zeros(length(phis),1);
    for i=1:length(theta)
        thetas(2*i^2-5*i+4:2*i^2-1*i) = theta(i);
    end
    E_exp = zeros(length(thetas),1);
    %for each sample location
    for i=1:length(thetas)
        E_exp(i) = E_expansion(Elmn,L,N,zeta,x(n),thetas(i),phis(i)); 
    end
    E(n,1:(L(n)+1)*(L(n)+2)*0.5) = E_exp;
end

