function Elmn_GT = generate_GT(Lmax,N,diff_time,zeta,REAL,varargin)
% Elmn_GT.m Find ground truth in SPF basis by using mutlishell scheme with 254 samples
%each shell and 10 shells - should be enough samples so good GT as lots
%more than degrees of freedom but will also have truncation error
%author: Alice Bates
% date last changed: 13/02/18

b_vec =[100,200,300,400,500,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250];
load('ESR_254.mat')
XYZB_big = [];
for i=1:length(b_vec)
    XYZB_temp = zeros(size(XYZB));
    XYZB_temp(:,1:3) = XYZB(:,1:3);
    XYZB_temp(:,4) = b_vec(i);
    XYZB_big = [XYZB_big;XYZB_temp];
end

if nargin == 7 %diffusivities single axisymmetric fibre
    [y,~,~] =  gen_Gauss_diff_sig(1e10,XYZB_big,varargin{1},varargin{2});%Gaussian mixture model with Rician noise added
elseif nargin == 6 %crossing angle
    [y,~,~] = gen_Gauss_diff_sig(1e10,XYZB_big,varargin{1});%Gaussian mixture model with Rician noise added
else
    [y,~,~] =gen_Gauss_diff_sig(1e10,XYZB_big);
end

%constructing SPF matrix
[basis] = construct_SPF_basis(Lmax, N, XYZB_big, diff_time,zeta,REAL);
S = transpose(basis);

Elmn_GT =  (ctranspose(S)*S)\(ctranspose(S))*y; 

