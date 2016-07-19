function [R,lag] = cross_corr_fun(tau, u, index)
%Computed the cross correlation as described in Marusic & Heuer 2007 PRL

test = false;
if test == false;
%t_long = [-t(end):-1:t(2); t];
loc = ceil(length(tau)/2);
m=1;
val = floor((length(tau)-loc)/index) - 1;
R = zeros(length(-index:1:index),1);
for i = -index:1:index;
    R(m) = mean(tau(loc-val*index:loc+val*index).*u(loc-val*index+i:loc+val*index+i))/...
        sqrt(mean(tau.^2)*mean(u.^2));
    m=m+1;
end
lag = [-index:1:index];
else
x=tau;
y=u;
R = max(abs(ifft(fft(x).*conj(fft(y)))))./...
    (sqrt(max(abs(ifft(fft(x).*conj(fft(x)))))).*...
    sqrt(max(abs(ifft(fft(y).*conj(fft(y)))))));
lag = [];
end


end

