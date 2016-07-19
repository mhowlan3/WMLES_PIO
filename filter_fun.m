function [ signalf ] = filter_fun(  signal , fcutoff , fs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
fs = fs;
fc_h = fcutoff;
fc_l = fcutoff;

filt_order  = 8;    % Set up Butterworth filter
flp  = fdesign.lowpass('N,F3dB', filt_order, fc_h, fs);
filter_lp = design(flp, 'butter'); % low pass filter
%fhp  = fdesign.highpass('N,F3dB', filt_order, fc_l, fs);
%filter_hp = design(fhp, 'butter'); % high pass filter

signalf = filtfilt(filter_lp.sosMatrix,filter_lp.ScaleValues,signal);

end

