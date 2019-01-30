function [] = nsht_plot_sampling(SHOW_ANTI)
% nsht_plot_sampling -  Plots the sampling points on the sphere over the
% surface of the sphere. For the spectral multi-shell sampling scheme.
% Plots first the actual samples and then the antipodal points in a lighter
% shade if SHOW_ANTI = 1, just plots actual samples otherwise
%author: Alice Bates
%date: 4/12/17

%parameters
LL=80;

[x, y, z] = sphere(LL);

figure('Color', [1 1 1]);
colourMap = colormap(gray);
brighten(gray,0.9);
h = surf(0.99*x,0.99*y,0.99*z,ones(size(x)));
set(h,'Linestyle', 'none');

hold on;
axis off;
axis equal;

%outer shell
L=9;
[THETA, FI] = nsht_sampling_points(L);


%%


%% longitudes

TT = [pi/5 2*pi/5 3*pi/5 4*pi/5];
FF = 0:0.01:2*pi;

for i=0:1:length(TT)-1
        Xm = 0.99*sin(TT(i+1)).*cos(FF);
        Ym = 0.99*sin(TT(i+1)).*sin(FF);
        Zm = 0.99*cos(TT(i+1)).*ones(size(FF));
        plot3(Xm,Ym,Zm, '--','Color',[0 0 0], 'linewidth',0.5);
end
%% latitudes
FF = pi/4:pi/4:2*pi;
TT = 0:0.01:pi;

for i=0:1:length(FF)-1
        Xm = 0.99*sin(TT).*cos(FF(i+1));
        Ym = 0.99*sin(TT).*sin(FF(i+1));
        Zm = 0.99*cos(TT);
        plot3(Xm,Ym,Zm, '--','Color',[0 0 0], 'linewidth',0.5);
        hold on;
end
% 



%%

for i=0:1:length(THETA)-1
    for j=2*i^2-i+1:2*i^2+3*i+1
        Xm = 1*sin(THETA(i+1)).*cos(FI(j));
        Ym = 1*sin(THETA(i+1)).*sin(FI(j));
        Zm = 1*cos(THETA(i+1));
       

        plot3(Xm,Ym,Zm, '.', 'markersize', 25,  'color',[0 0 1]);
    
        if SHOW_ANTI
            Xm_anti = -Xm;
            Ym_anti = -Ym;
            Zm_anti = -Zm;
            plot3(Xm_anti,Ym_anti,Zm_anti, '.', 'markersize', 25,  'color',[.7 .7 1]);
        end

    end
end



%% second largest shell
L=7;
[THETA, FI] = nsht_sampling_points(L);
%%


%%

for i=0:1:length(THETA)-1
     for j=2*i^2-i+1:2*i^2+3*i+1
        Xm = 1*sin(THETA(i+1)).*cos(FI(j));
        Ym = 1*sin(THETA(i+1)).*sin(FI(j));
        Zm = 1*cos(THETA(i+1));
   
        plot3(Xm,Ym,Zm, '.', 'markersize', 25,  'color',[1 0 0]);
        if SHOW_ANTI
            Xm_anti = -Xm;
            Ym_anti = -Ym;
            Zm_anti = -Zm;
            plot3(Xm_anti,Ym_anti,Zm_anti, '.', 'markersize', 25,  'color',[1 0.8 0.8]);
        end
     end
end


%% second smallest shell
L=5;
[THETA, FI] = nsht_sampling_points(L);

%%

for i=0:1:length(THETA)-1
    for j=2*i^2-i+1:2*i^2+3*i+1
        Xm = 1*sin(THETA(i+1)).*cos(FI(j));
        Ym = 1*sin(THETA(i+1)).*sin(FI(j));
        Zm = 1*cos(THETA(i+1));
        plot3(Xm,Ym,Zm, '.', 'markersize', 25,  'color','g');

        if SHOW_ANTI
            Xm_anti = -Xm;
            Ym_anti = -Ym;
            Zm_anti = -Zm;
            plot3(Xm_anti,Ym_anti,Zm_anti, '.', 'markersize', 25,  'color',[0.9 1 0.9]);
        end

    end
end



%% smallest shell
L=3;
[THETA, FI] = nsht_sampling_points(L);

for i=0:1:length(THETA)-1
     for j=2*i^2-i+1:2*i^2+3*i+1
        Xm = 1*sin(THETA(i+1)).*cos(FI(j));
        Ym = 1*sin(THETA(i+1)).*sin(FI(j));
        Zm = 1*cos(THETA(i+1));
   
        plot3(Xm,Ym,Zm, '.', 'markersize', 25,  'color','k');
        
        if SHOW_ANTI
            Xm_anti = -Xm;
            Ym_anti = -Ym;
            Zm_anti = -Zm;
            plot3(Xm_anti,Ym_anti,Zm_anti, '.', 'markersize', 25,  'color', colourMap(40,:));
        end

    end
end


end


