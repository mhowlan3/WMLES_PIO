function [] = PIO_cali_fast(d, plane_num,retau, DNS_data, save_vals)
%Calibration file for the Predictive Inner-Outer model of Marusic et al.
%2011 with DNS data


%% Unpack the DNS data
planes_tau_xy = DNS_data.planes_tau_xy;
planes_U = DNS_data.planes_U;
load(strcat('sim_param', num2str(retau),'.mat'))

%mkdir
if exist(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/'))~=7;
    mkdir(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/'));
end

disp('Retau='); disp(Retau);
disp('Filter length scale in delta ='); disp(d);
disp('Velocity plane number='); disp(plane_num);

%% Preprocess
u_filt = zeros(size(planes_tau_xy));
t_filt = zeros(size(planes_tau_xy));
tau_wall_filt = zeros(size(planes_tau_xy));
tau_t_filt = zeros(size(planes_tau_xy));
fplus_tau = zeros(size(planes_tau_xy(1,:,:)));
fplus_u = zeros(size(planes_tau_xy(1,:,:)));
tau_wall_unfilt = zeros(size(planes_tau_xy));
disp('Start calibration data collection')
nx = length(planes_tau_xy(:,1,1));
nz = length(planes_tau_xy(1,:,1));
nt = length(planes_tau_xy(1,1,:));
m_store = [];
xm_store = [];
theta_store = [];
alpha_store = [];
rms_val_store = [];
load(strcat('sim_param', num2str(retau),'.mat'))
for t = 1 : nt;
    %%%Filter length scale for low pass filter keeping only large scales
    %U_c = mean(mean(planes_U(:,:,plane_num,t)));
    f_filt = (U_c(plane_num, t)/utau)/(d*Retau);
    %f_filt_tau = ((12.1*utau)/utau)/(d*Retau);
    f_filt_tau=f_filt;
     lambda_x = (d*Retau*nu)/utau;
    %%%
    for i = 1 : nz;
        %tic
    %Wall shear stress
    [series_tau_xy, time_tau_xy, ~] = ...
        taylorshyp(planes_tau_xy(:,i,t), utau, x, false, U_c(plane_num, t), lambda_x);
    series_tau_xy = series_tau_xy - mean(series_tau_xy);
    tplus = time_tau_xy*(utau^2)/nu;
    fplus_tau(i,t) = 1/mean(diff(tplus));
    cutt = filter_fun(series_tau_xy', f_filt_tau, fplus_tau(i,t)); %2.65*10^-3
    
    tau_wall_filt(:,i,t) = cutt/(utau^2);
    tau_t_filt(:,i,t) = tplus;
    tau_wall_unfilt(:,i,t) = series_tau_xy/(utau^2);
    
    %Velocity
    u_space = planes_U(:,i,plane_num,t);
    [u_series, t_series, ~] = taylorshyp(u_space, utau, x, false, U_c(plane_num,t), lambda_x);
    u_series = u_series - mean(u_series);
    %figure
    tplus = t_series*(utau^2)/nu;
    fplus_u(i,t) = 1/mean(diff(tplus));
    cut_u = filter_fun(u_series', f_filt, fplus_u(i,t));
    u_filt(:,i,t) = cut_u/utau;
    t_filt(:,i,t) = tplus;
    
    %%%% 
    x_corr = [2, 1]; %x limits in delt for the cross correlation computation
    % x_corr = 2 seems to give good values for xm and theta

    t_cutoff_vect = find(tau_t_filt(:,i,t) > t_filt(end,i,t));
     if numel(t_cutoff_vect)==0;
         t_cut_ind = length(tau_t_filt(:,i,t));
     else
        t_cut_ind = t_cutoff_vect(1)-1; %-1
     end
     t_corr = tau_t_filt(1:t_cut_ind,i,t);
     u_corr = interp1(t_filt(:,i,t), u_filt(:,i,t), t_corr);   
     
     t_val = x_corr./U_c(plane_num,t);
     t_val_plus = t_val*(utau^2)/nu;
     tmp = abs(t_corr-t_val_plus(1));
     [idx index(1)] = min(tmp);
     tmp = abs(t_corr-t_val_plus(2));
     [idx index(2)] = min(tmp);
     %disp(index(1))  
     %tic
%      [cross_corr{i,t},lag{i,t}] = cross_corr_fun(tau_wall_filt(1:t_cut_ind,i,t), u_corr,index(1));
     [cross_corr{i,t}, lag{i,t}] = xcorr(u_corr, tau_wall_filt(1:t_cut_ind,i,t), 'coeff');
     ind1 = find(lag{i,t} == -index(1));
     ind2 = find(lag{i,t} == index(1));
     cross_corr{i,t} = cross_corr{i,t}(ind1:ind2);
     lag{i,t} = lag{i,t}(ind1:ind2);
%      figure
%      hold on
%      plot(lag{i,t}, cross_corr{i,t})
%      plot(lag2{i,t}, cross_corr2{i,t})
%      [max_val, loc] = max(cross_corr{i,t})
%      [max_val2, loc2] = max(cross_corr2{i,t})
     %toc
     %[cross_corr{s,t},lag{s,t}] = cross_corr_fun(tau_wall_filt(:,s,t), u_filt(:,s,t),index(1));
    hold on
    %plot(lag{s,t}, cross_corr{s,t})
     
     [m,max_loc] = max(cross_corr{i,t});
     alpha = m.*(sqrt(mean(tau_wall_filt(1:t_cut_ind,i,t).^2))./sqrt(mean(u_corr.^2)));
     rms_val = (sqrt(mean(tau_wall_filt(1:t_cut_ind,i,t).^2))./sqrt(mean(u_corr.^2)));
    x_long = [-x(t_cut_ind:-1:2);x(1:t_cut_ind)];
    if lag{i,t}(max_loc)==0;
        lag{i,t}(max_loc)=1;
    end   
     tt = t_corr(abs(lag{i,t}(max_loc)));
     tt = tt*fplus_tau(i,t)/fplus_u(i,t); %convert into the log layer time scale
     xx = tt*(U_c(plane_num,t)/utau);
     xm = xx*(nu/utau);
     %xm = x(abs(lag{s,t}(max_loc)));
     %xm= x_long(max_loc);
     theta = atand(ypos(plane_num)/(xm));
    m_store = [m_store; m];
    alpha_store = [alpha_store; alpha];
    xm_store = [xm_store; xm];
    theta_store = [theta_store; theta];
    rms_val_store(i) = rms_val;
    %toc
    end
    rms_mean_val(t) = mean(rms_val_store);
    clear rms_val_store
    disp(t)
end
if save_vals == true;
save(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/', 'cali_vars.mat'), 'u_filt', 't_filt', 'tau_wall_filt',...
    'tau_wall_unfilt', 'tau_t_filt', 'fplus_tau', 'fplus_u', '-v7.3')
save(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/', 'tau_wall_unfilt.mat'), 'tau_wall_unfilt', 'tau_t_filt', 'fplus_tau', 'fplus_u', '-v7.3')
end
disp('Preprocess finished')

%% Calibrate
disp('Start calibration')
ss=size(cross_corr);
xm_mean = zeros(nt,1);
alpha_mean = zeros(nt,1);
theta_mean = zeros(nt,1);
for k = 1 : nt;
    corr_sum = zeros(size(cross_corr{1,k}));
    for i = 1 : nz;
        corr_sum = corr_sum+cross_corr{i,k};
    end
    cross_corr_mean{k} = corr_sum/ss(1);
    [m,max_loc] = max(cross_corr_mean{k});
    %xm_mean(k) = x(abs(lag{1,k}(max_loc))); % if using spatial correlation
    tt = t_corr(abs(lag{i,t}(max_loc)));
    tt = tt*fplus_tau(i,t)/fplus_u(i,t); %convert into the log layer time scale
    xx = tt*(U_c(plane_num,t)/utau);
    xm_mean(k) = xx*(nu/utau);
    alpha_mean(k) = m.*rms_mean_val(k);
    theta_mean(k) = atand(ypos(plane_num)/xm_mean(k));
end    
theta_avg = atand(ypos(plane_num)./(xm_mean));


%Mean values from the constants calculated at each spanwise location
alpha_std = std(alpha_store);
theta_std = std(theta_store);
% 
% %Save data
if save_vals == true;
save(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
   num2str(d), 'd/', 'cali_avg_vals.mat'), 'alpha_store', 'theta_store',...
   'alpha_mean', 'alpha_std', 'theta_mean', 'theta_std', 'xm_mean', 'theta_avg',...
   'xm_store', 'cross_corr_mean', 'cross_corr', 'lag', '-v7.3')
end
% disp('Cali done')
% disp(strcat('plane_',num2str(plane_num),'/', num2str(d), 'd'))
disp('Mean calibration done')

%% Universal Velocity signal

% Reconstruct the unfiltered tau signal to correspond to the lower
% convection velocity at that location
load(strcat('sim_param', num2str(retau),'.mat'))
for t =  1 : nt;
    tau_t_filt(:,:,t) = tau_t_filt(:,:,t) * U_c(plane_num,t)/(12.1*utau);
    fplus_tau(1:nz,t) = 1/mean(diff(tau_t_filt(:,1,t)));
end
clear tau_w_star
disp('Start universal signal calibration')
for t = 1 : nt;
    for i = 1 : nz;
        %Unpack some calibration values
%         alpha(i,t) = ...
%             alpha_store(i+((t-1)*length(tau_wall_unfilt(1,1,:)))); %this uses the xm shift per spanwise location
        %Use mean xm or individual spanwise location
        xm(i,t) = xm_mean(t); %%% this uses the xm shift per plane snapshot
        alpha(i,t) = alpha_mean(t); % this uses the alpha per plane snapshot
        %xm(i,t) = mean(xm_mean(:));
        %xm(i,t) = xm_store(i+((t-1)*length(tau_wall_unfilt(1,1,:)))); %this uses the xm shift per spanwise location
        if xm(i,t) == 0
            xm(i,t) = x(2);
        end
           
        delta_t_avg = mean(diff(t_filt(:,i,t)));
        shift(i,t) = ceil((((abs(xm(i,t))/U_c(plane_num,t))*(utau^2)/nu)/delta_t_avg));
        %shift_round(i,t) = ceil(shift(i,t));
        
        %%% u_OL, the large scale motion velocity in the outer log layer
        %%% must be shifted forward in the streamwise direction to match
        %%% with the shear stress signal at the wall
        %disp(shift(i,t))
        u_OL{i,t} = u_filt(shift(i,t):nx,i,t);
        t_cutoff_vect = find(tau_t_filt(:,i,t) > t_filt(end,i,t)-((shift(i,t)+1)*delta_t_avg));
        if numel(t_cutoff_vect)==0;
            t_cut_ind(i,t) = length(tau_t_filt(:,i,t));
        else
            t_cutoff(i,t) = t_cutoff_vect(1); %-1
        end
        %t_cutoff(i,t) = t_cutoff_vect(1);
        u_OL_interp{i,t} = interp1(fplus_u(i,t)*t_filt(1:end-shift(i,t)+1,i,t)/fplus_tau(i,t), u_OL{i,t}, tau_t_filt(1:t_cutoff(i,t)));
        % the multiplication of t_filt by (fplus_u/fplus_tau) shifts the
        % u_OL parameter from the temporal space of the log layer to that
        % of the wall
        t_univ{i,t} = tau_t_filt(1:t_cutoff(i,t));
        
        %Universal wall shear stress
        tau_w_star{i,t} = (tau_wall_unfilt(1:t_cutoff(i,t),i,t)...
            - alpha(i,t)*u_OL_interp{i,t}')./...
            (1 + alpha(i,t)*u_OL_interp{i,t}');
        
        %any(isnan(tau_w_star{i,t}))
    end
end
disp('Universal signal calibration done')

%% Tau specta
fplus = 1./mean(diff(tau_t_filt(:,1,1)));
ss=size(tau_w_star);
nw=5;
disp('start tau_w_star loop')
for t = 1 : nt;
    for j = 1 : nz;
%         xx=tau_w_star{j,t}'; ns=floor(length(tau_w_star{j,t})/nw);
%         [pxx(:,j,t), ff(:,j,t), ~] = pwelch(xx,hanning(length(tau_w_star{j,t})),[],[],mean(fplus));
         kk = length(tau_w_star{j,t});
%         [pxx2(:,j,t), ff2(:,j,t), ~] = pwelch(tau_wall_unfilt(1:kk,j,t),hanning(length(tau_w_star{j,t})),[],[],mean(fplus));
        Fs = fplus;
        x = tau_w_star{j,t}';
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:N/2+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/length(x):Fs/2;
        if j > 1 || t>1;
        if length(freq)==length(ff(:,1,1));
        ff(:,j,t)=freq;
        pxx(:,j,t)=psdx;
        else
            ff(:,j,t) = ff(:,1,1);
            pxx(:,j,t) = interp1(freq, psdx, ff(:,1,1));
            pxx(end,j,t) = 0;
        end
        else
            ff(:,j,t)=freq;
            pxx(:,j,t)=psdx;
        end
        
        %%%
        x = tau_wall_unfilt(1:kk,j,t);
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:N/2+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/length(x):Fs/2;
        if j > 1 || t>1;
        if length(freq)==length(ff2(:,1,1));
        ff2(:,j,t)=freq';
        pxx2(:,j,t)=psdx';
        else
            ff2(:,j,t) = ff2(:,1,1);
            pxx2(:,j,t) = interp1(freq', psdx', ff2(:,1,1));
            pxx2(end,j,t) = 0;
        end
        else
            ff2(:,j,t)=freq;
            pxx2(:,j,t)=psdx;
        end
        %disp(j)
    end
    disp(t)
end
%this test is using pwelch with no window averaging
%save('test_tau_w_spectra_', 'pxx', 'ff','pxx2','ff2', '-v7.3')
%end
f1=figure;
hold on
px = mean(pxx,3);
px = mean(px, 2);
px2 = mean(pxx2,3);
px2 = mean(px2, 2);
f = mean(ff, 3);
f = mean(f,2);
semilogx(f, px.*f);
semilogx(f, px2.*f);
set(gca, 'xscale', 'log')
%set(L, 'interpreter', 'latex')
xlabel('$f^{+}$', 'interpreter', 'latex','fontsize', 20)
ylabel('$f^{+} \phi_{\tau^{+} \tau^{+}}$', 'interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20)
L=legend('$\tau_{w}^{*}$','$\tau_{w}$');
set(L,'interpreter','latex');
box on
title(strcat('$ x/\delta = ', num2str(d), '$'), 'interpreter', 'latex', 'fontsize', 20) 

%% Save the final results
disp('Start save')
if exist(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/'))~=7;
    mkdir(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/'));
end
save(strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/', 'cali_param_final.mat'),...
    'alpha_mean', 'theta_mean', 'tau_w_star', 'tau_wall_unfilt',...
    't_univ', 'f', 'px', 'ff', 'pxx', 'pxx2', '-v7.3')
disp('Save finished')
%save power spectra figure
saveas(f1, strcat('data/', 'Retau',num2str(retau),'/','plane_',num2str(plane_num),'/',...
    num2str(d), 'd/', 'tau_w_star_spectra'), 'fig')

end