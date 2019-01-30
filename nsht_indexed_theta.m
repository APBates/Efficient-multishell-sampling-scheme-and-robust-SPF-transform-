function thetas = nsht_indexed_theta(L)
%place (L+1)/2 equiangular rings around theta. Order these so that
%descending order from closest to equator

TT_temp = pi*(2*(0:(L-1)/2)+1)/L;

[ ~, T_index] = sort(abs(TT_temp-pi/2),'descend');
thetas = TT_temp((T_index));


