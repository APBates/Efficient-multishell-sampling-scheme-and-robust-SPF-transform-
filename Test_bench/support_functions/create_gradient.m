
function [] = createGradient(bvalue,L)
%create_gradient.m used to create the Nx4 sampling matrix where N = number
%of sample points each row (x,y,z,b-value)
B_VALUE = bvalue/1e6;%411.3;
NUM_SAMPLES = (L+1)*L/2;
[theta, phi] = nsht_sampling_points(L);
XYZB = zeros(NUM_SAMPLES, 4);

XYZB(:,4) = B_VALUE;
theta_expanded = zeros(size(phi));
for i = 1:length(theta)
    theta_expanded(2*i^2-5*i+4:2*i^2-i)=theta(i);
end

for i = 1: length(phi)
    XYZB(i,1) = sin(theta_expanded(i))*cos(phi(i));
    XYZB(i,2) = sin(theta_expanded(i))*sin(phi(i));
    XYZB(i,3) = cos(theta_expanded(i));
end

save(['gradient_optimal_L' num2str(L) '_B' num2str(B_VALUE) '.mat'],'XYZB');



