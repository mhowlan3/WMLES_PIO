function [u_time, time, U_c] = taylorshyp(u_space, utau, x, tau, U_c, lambda_x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
len = length(u_space);
x_diff = diff(x);

%vectorize
if abs(mean(diff(x_diff))) < 1e-20;
    if tau == true;
        u_time = u_space(end:-1:1)';
        start_t = 0;
        end_t = (x(end)-x(1))/(12.1*utau);
        time = linspace(start_t, end_t, len);
        U_c = NaN;
    else
        u_time = u_space(end:-1:1)';
        start_t = 0;
        end_t = (x(end)-x(1))/(U_c);
        time = linspace(start_t, end_t, len);
    end
    
    
else
if tau == true;
   
f = 12.1*(utau/lambda_x);

for i = 1 : len;
    u_time(len+1-i) = u_space(i);
    time(len+1-i) = (x(end)-x(i))/(12.1*utau);
    %time(len+1-i) = (x(end)-x(i))/(U_c);
end
U_c = NaN;
else
    %U_c = mean(u_space);
    for i = 1 : len;
        u_time(len+1-i) = u_space(i);
        time(len+1-i) = (x(end)-x(i))/(U_c);
        
        %time(len+1-i) = (i*7000)/(12.1*utau);
    end
    f = 1./diff(time);
end
end


end

