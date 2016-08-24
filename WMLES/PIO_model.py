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


# Input parameters
d = 2.15
retau = 2000
plane_num = 1
folder = ''
folder = '/home/mhowland/PIO_model/data/' + 'Retau' + str(retau) + '/' + 'plane_' + str(plane_num) + '/' + str(d) + 'd/'


# Some operation functions
def rms(val):
        rms_val = np.sqrt(np.mean(np.square(val)))
        return rms_val


# Load calibration data
nt = 10
for t in range(1, nt+1):
        nz = 4608 #number of spanwise locations
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
#print(tau_w_star.shape)
print('Calibration data loaded')

#print(min_len)
# Import the WMLES velocity data
folder_bin = 'WMLES_u_binary.bin'
f = open(folder_bin, 'rb')
data = np.fromfile(f, dtype=np.float64)
#print(len(data))
#WMLES dataset size
# sorry this is currently hardcoded but just take the size from the WMLES file 
u_size = (11400, 250, 160)
u_size = (11400, 2, 2)
u_WMLES = np.zeros(u_size, dtype=np.float64)
for k in range(0, u_size[2]):
	for j in range(0,u_size[1]):
		for i in range(0, u_size[0]):
			u_WMLES[i,j,k] = data[i+ j*u_size[1] + k*u_size[2]] 
folder_bin = 'WMLES_t_binary.bin'
f = open(folder_bin, 'rb')
t_WMLES = np.fromfile(f, dtype=np.float64)
#print(t_WMLES.shape) 
print('WMLES data loaded')
#t_WMLES = t_WMLES[0:1000]
#print(t_WMLES.shape)
#print(u_WMLES.shape) 

# Apply PIO model function
def PIO_model_apply(u_WMLES_vect, tau_w_star_vect, t_WMLES_vect, t_tau_vect, d, retau, alpha):
	f_wmles = 1.0/np.mean(np.diff(t_tau_vect))
	U_c_plus = np.mean(u_WMLES_vect)
	u_OL = u_WMLES_vect - U_c_plus
	#print(t_WMLES_vect.shape)
	#print(u_WMLES_vect.shape)
	u_OL = np.interp(t_tau_vect, t_WMLES_vect, u_WMLES_vect)	
	f_filt = (U_c_plus)/(float(d*retau))
	u_OL = butter_lowpass_filter(u_OL, f_filt, f_wmles, 8)
	tau_w_predic = np.zeros(len(tau_w_star_vect), dtype=np.float64) 
	for i in range(0, len(tau_w_star_vect)-1):
		tau_w_predic[i] = tau_w_star_vect[i]*(1+alpha*u_OL[i]) + alpha*u_OL[i]
	return tau_w_predic, f_wmles

#Filtration definition, lowpass butterworth filter
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

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
cali_size = tau_w_star.shape
#print(cali_size)
#cali_size = [len(tau_w_star[:,0,0]), 3, 2]
tau_wp = np.zeros(cali_size, dtype=np.float64)
folder_bin = folder + 'alpha.bin'
f = open(folder_bin, 'rb')
alpha = np.fromfile(f, dtype=np.float64)
N = len(tau_w_star[:,0,0])
psdx_size = [float(N)/2.0+1.0, cali_size[1], cali_size[2]]
tau_wp_psdx = np.zeros(psdx_size, dtype=np.float64)
tau_wp_f = np.zeros(psdx_size, dtype=np.float64)
print(alpha[0], alpha[1])

for j in range(0, cali_size[1]):
	for k in range(0, cali_size[2]):
		tau_wp[:,j,k], Fs = PIO_model_apply(u_WMLES[:,ii,jj], tau_w_star[:,j,k],t_WMLES, t_tau[:, k], d, retau, alpha[k])
		#print(tau_wp_psdx.shape)
		#wut1,wut2 = psdx(tau_wp[:,j,k],Fs)
		#print(wut1.shape)
		tau_wp_psdx[:,j,k], tau_wp_f[:,j,k] = psdx(tau_wp[:,j,k], Fs)	 
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
#print(tau_wp_psdx_mean[0])
#print(tau_wp_psdx_mean[100])
# Plots
print('Begin plotting')
#from matplotlib import rc

#
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.unicode'] = True
#


#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
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
plt_name = folder + 'spectra_fig.eps'
plt.savefig(plt_name, format='eps', dpi=1000)
print('Plots saved')
#plt.show()

# Output
sys.stderr.close()
sys.stderr = sys.__stderr__




