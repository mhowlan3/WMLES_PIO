#!/usr/bin/env python

##############################################
### Implement PIO model with WMLES results ###
### Michael Howland, mhowland@stanford.edu ###
### 8/16/16                                ###
##############################################

## Import some libraries
import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib as mpl
mpl.use('Agg')
#mpl.use('GTK')
#mport matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import sys
sys.stderr = open('errorlog.txt', 'w')
from matplotlib import rc

# Input parameters
d = 2.5
retau = 2000
plane_num = 1
folder = ''
#WM_type = 'integral'
WM_type = 'EQWM'
folder = '/home/mhowland/PIO_model/data/' + 'Retau' + str(retau) + '/' + 'plane_' + str(plane_num) + '/' + str(d) + 'd/'


# Some operation functions
def rms(val):
        rms_val = np.sqrt(np.mean(np.square(val)))
        return rms_val

# Load calibration data
# Load the size of the calibration data
if retau == 2000:
	nt = 10
	nz = 4608
else:
	nt = 20
	nz = 3072 
for t in range(1, nt+1):
        #nz = 4608 #number of spanwise locations
        folder_bin = folder + 'cali_mat' + '_t' + str(t) + '.bin'
        #print folder_bin
        f = open(folder_bin, 'rb')
        data = np.fromfile(f, dtype=np.float64)
        #print len(data)
        nx = float(len(data))/(float(nz)+1)
        #print nx
	if t == 1:
        	s = (nx, nt)
        	tplus = np.zeros(s, dtype=np.float64)
        	s = (nx, nz, nt)
		ss = (nz, nt)
        	tau_w_star = np.zeros(s, dtype=np.float64)
        	tau_rms = np.zeros(ss, dtype=np.float64)
		min_len = nx
	if nx < min_len:
		min_len = nx
	# Since the DNS data has varying convection velocities it is possible
	# that the time series may carry differing lengths, to account for this
	# only use up to the value of min_len stored here when applying 
	# tau_w_star
        # Unpack the binary data
        for j in range(0, nz):
		#print(j)
                if j == 0:
			temp = data[0 : nx]
			#print(nx)
			#print(len(temp))
			if len(temp) < s[0]:
				for c in range(0, int(s[0]-len(temp))):
					temp = np.append(temp, 0.0)
			tplus[:,t-1] = temp[0 : s[0]]
                else:
			temp = data[(j)*nx : (j+1)*nx]
                        if len(temp) < s[0]: 
                                for c in range(0, int(s[0]-len(temp))):
                                        temp = np.append(temp, 0.0)   
			tau_w_star[:,j-1,t-1] = temp[0 : s[0]]
			temp = rms(data[(j)*nx : (j+1)*nx]) 
                        tau_rms[j-1,t-1] = temp
tau_w_star = tau_w_star[0:min_len, :, :]
t_tau = tplus[0:min_len, :]
print('Calibration data loaded')

# Import the WMLES velocity data
def read_bin(file):
	f = open(file, 'rb')
	bin_data = np.fromfile(f, dtype=np.float64)
	return bin_data
if WM_type == 'integral':
	folder_bin = 'integral_WMLES_time.bin'	
	t_WMLES_1 = read_bin(folder_bin)
	t_WMLES = t_WMLES_1 - t_WMLES_1[0] # normalize b/c Xiang's time didn't start from 0
	nu_xiang = 1.0/2000.0
	t_WMLES = [xxx * 2000.0 for xxx in t_WMLES] # nondimensionalize in inner variables
	folder_bin = 'integral_WMLES_u_series.bin'
	u_WMLES = read_bin(folder_bin)
	height = '1'
else: 
	height = '5'
	folder_bin = 'WMLES_u_binary_y' + height + '.bin'
	data = read_bin(folder_bin)
	# WMLES dataset size
	# sorry this is currently hardcoded but just take the size from the WMLES file 
	u_size = (11400, 250, 160)
	u_size = (11400, 2, 2)
	u_WMLES = np.zeros(u_size, dtype=np.float64)
	for k in range(0, u_size[2]):
		for j in range(0,u_size[1]):
			for i in range(0, u_size[0]):
				u_WMLES[i,j,k] = data[i+ j*u_size[1] + k*u_size[2]]
	folder_bin = 'WMLES_t_binary_y' + height + '.bin'
	t_WMLES = read_bin(folder_bin)
print('WMLES data loaded')


# Filtration definition, lowpass butterworth filter
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

# Apply PIO model function
def PIO_model_apply(u_WMLES_vect, tau_w_star_vect, t_WMLES_vect, t_tau_vect, d, retau, alpha):
	f_wmles = 1.0/np.mean(np.diff(t_tau_vect))
	U_c_plus = np.mean(u_WMLES_vect)
	u_OL_vect = u_WMLES_vect - U_c_plus
	#print(t_WMLES_vect.shape)
	#print(u_WMLES_vect.shape)
	u_OL_temp = np.interp(t_tau_vect, t_WMLES_vect, u_OL_vect)	
	f_filt = (U_c_plus)/(float(d*retau))
	u_OL = butter_lowpass_filter(u_OL_temp, f_filt, f_wmles)
	tau_w_predic = np.zeros(len(tau_w_star_vect), dtype=np.float64) 
#	## test for debugging spectra problem
#	fig3 = plt.gca()
#	fig3.plot(t_tau_vect, tau_w_star_vect, 'b-', t_tau_vect, u_OL_temp, 'k-', t_tau_vect, u_OL, 'g-')
#	plt_name_3 = folder + 'debug_plot.eps'
#	plt.savefig(plt_name_3, format='eps', dpi=1000)
#	plt.hold(False)
#	print(f_filt)
#	##
	for i in range(0, len(tau_w_star_vect)-1):
		tau_w_predic[i] = tau_w_star_vect[i]*(1+alpha*u_OL[i]) + alpha*u_OL[i]
	return tau_w_predic, f_wmles, u_OL

# Power spectral density
def psdx(x, Fs):
	N = len(x)
	xdft = np.fft.fft(x)
	xdft = xdft[0:float(N)/2+1]
	psdx = (1/(Fs*N)) * np.square(abs(xdft))
	psdx[1:len(psdx)-1] = 2*psdx[1:len(psdx)-1]
	#print(Fs)
	#print(Fs/N)
	#print(N)
	freq = np.arange(0, Fs/2, Fs/(N+1))
	return psdx, freq

# apply the PIO model to the data WMLES dataset
ii=1
jj=1
## ii and jj dictate the spanwise and time locations of the WMLES dataset, results may vary
cali_size = tau_w_star.shape
#print(cali_size)
#cali_size = [len(tau_w_star[:,0,0]), 3, 2]
tau_wp = np.zeros(cali_size, dtype=np.float64)
u_OL_mat = np.zeros(cali_size, dtype=np.float64)
folder_bin = folder + 'alpha.bin'
f = open(folder_bin, 'rb')
alpha = np.fromfile(f, dtype=np.float64)
N = len(tau_w_star[:,0,0])
psdx_size = [float(N)/2.0+1.0, cali_size[1], cali_size[2]]
tau_wp_psdx = np.zeros(psdx_size, dtype=np.float64)
tau_wp_f = np.zeros(psdx_size, dtype=np.float64)
u_OL_psdx = np.zeros(psdx_size, dtype=np.float64)
u_OL_f = np.zeros(psdx_size, dtype=np.float64)
print(alpha[0], alpha[1])

for j in range(0, cali_size[1]):
	for k in range(0, cali_size[2]):
		if WM_type == 'integral':
			#print('this if statement is working')
			u_WMLES_apply = u_WMLES
		else:
			#print('this if statement is NOT working')
			#print(j, k)
			u_WMLES_apply = u_WMLES[:,ii,jj]
		tau_wp[:,j,k], Fs, u_OL_mat[:,j,k] = PIO_model_apply(u_WMLES_apply, tau_w_star[:,j,k],t_WMLES, t_tau[:, k], d, retau, alpha[k])
		#print(tau_wp_psdx.shape)
		#wut1,wut2 = psdx(tau_wp[:,j,k],Fs)
		#print(wut1.shape)
		tau_wp_psdx[:,j,k], tau_wp_f[:,j,k] = psdx(tau_wp[:,j,k], Fs)
		u_OL_psdx[:,j,k], u_OL_f[:,j,k] = psdx(u_OL_mat[:,j,k], Fs)	 
###########
# Make a plot to see the difference in the scales of fluctuation between u and tau
#print(rms(alpha[0]*u_WMLES[:,ii,jj]))
#print(rms(u_WMLES[:,ii,jj]))
#print(rms(tau_w_star[:,0,0]))

fig2 = plt.gca()
if WM_type == 'integral':
	u_WMLES_plot = u_WMLES
else:
	u_WMLES_plot = u_WMLES[:,ii,jj]
u_WMLES_plot = u_WMLES_plot - np.mean(u_WMLES_plot)
fig2.plot(t_WMLES, alpha[0]*u_WMLES_plot, 'b-', t_tau[:, 0], tau_w_star[:,0,0], 'k-')
plt_name_2 = folder + 'show_mag_diff' + WM_type + '.eps'
plt.savefig(plt_name_2, format='eps', dpi=1000)
plt.hold(False)

print(rms(tau_wp[:,0,0]))
#print(tau_wp.shape)
print('Model applied and psdx computed')

#Mean PSDX
#print(tau_wp_psdx.shape)
#print(tau_wp_f[0,0,0])


summ = np.zeros(psdx_size[0], dtype=np.float64)
for j in range(0, psdx_size[1]):
	for k in range(0, psdx_size[2]):
		summ = summ + tau_wp_psdx[:,j,k]
tau_wp_psdx_mean = summ/(float(j*k))


# Plots
print('Begin plotting')
fig1 = plt.gca()
psdx_plot = np.zeros(tau_wp_psdx_mean.shape, dtype=np.float64)
for i in range(0, len(tau_wp_psdx_mean)):
	psdx_plot[i] = (tau_wp_f[i,0,0] * tau_wp_psdx_mean[i])/1.035**2
# Plot the DNS spectra
f = open('DNS_tau_spectra.bin', 'rb')
DNS_tau_spectra = np.fromfile(f, dtype=np.float64)
DNS_f = DNS_tau_spectra[0:22]
DNS_psdx = DNS_tau_spectra[22:44]
#print(psdx_plot.shape)
#print(tau_wp_f[0,0,0])
#print(tau_wp_f[100,0,0])
#print(psdx_plot[100])
fig1.plot(2*np.pi*tau_wp_f[:,0,0], psdx_plot, 'b-', DNS_f, DNS_psdx, 'k-')

fig1.set_yscale('log')
fig1.set_xscale('log')
plt.ylim([10**-4, 10**0])
plt.xlim([10**-3, 2**1])
plt.xlabel(r"$w^{+}$", fontsize=16)
plt.ylabel(r"$w^{+} E(w^{+})/ \tau_{w}^{+2}$",fontsize=16)
plt.legend([r"$\tau_{wp, WMLES}^{\prime +}$", r"$\tau_{DNS}^{\prime +}$"], 'upper left')
title_str = r'$ \lambda_{x}^{+}=' + str(retau*d) + '$'
plt.title(title_str)
plt_name = folder + 'spectra_fig_y' + height + '_'  +  WM_type + '.eps'
plt.savefig(plt_name, format='eps', dpi=1000)

### Make u_OL PSDX plot ###
fig4 = plt.gca()
summ = np.zeros(psdx_size[0], dtype=np.float64)
for j in range(0, psdx_size[1]):
        for k in range(0, psdx_size[2]):
                summ = summ + u_OL_psdx[:,j,k]
u_OL_psdx_mean = summ/(float(j*k))
psdx_plot = np.zeros(u_OL_psdx_mean.shape, dtype=np.float64)
for i in range(0, len(u_OL_psdx_mean)):
        psdx_plot[i] = u_OL_psdx_mean[i]
fig4.plot(2*np.pi*u_OL_f[:,0,0], psdx_plot, 'b-')
fig4.set_yscale('log')
fig4.set_xscale('log')
#plt.ylim([10**-4, 10**0])
#plt.xlim([10**-3, 2**1])
plt.xlabel(r"$w^{+}$", fontsize=16)
plt.ylabel(r"$E(w^{+})$",fontsize=16)
title_str = r'$ \lambda_{x}^{+}=' + str(retau*d) + '$'
plt.title(title_str)
plt_name = folder + 'u_OL_spectra_' + height + '_'  +  WM_type + '.eps'
plt.savefig(plt_name, format='eps', dpi=1000)
###########################

print('Plots saved')
#plt.show()

#Write rms value to text file
#import os.path
rms_mean_sum = 0.0
for j in range(0, cali_size[1]):
        for k in range(0, cali_size[2]):
		rms_mean_sum = rms_mean_sum + rms(tau_wp[:,j,k])
rms_mean = rms_mean_sum / (int(cali_size[1]*cali_size[2]))
file_rms = "rms_data.txt"
# Appends the file, assumes that the file exists already
target = open('rms_data.txt', 'a')
target.write('%s %s %s\n' % (retau, d, WM_type))
target.write('Tau_wp_rms =  %s\n' % (rms_mean))
target.close()



# Output
sys.stderr.close()
sys.stderr = sys.__stderr__




