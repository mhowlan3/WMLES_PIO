#!/usr/bin/env python

##############################################
### Implement PIO model with WMLES results ###
### Michael Howland, mhowland@stanford.edu ###
### 8/16/16                                ###
##############################################


## Import some libraries
import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt
import sys
sys.stderr = open('errorlog.txt', 'w')


# Input parameters
d = 3.5
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
u_size = (11400, 250, 160)
#u_size = (11400, 2, 2)
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
	xdft = xdft[1:float(N)/2+1]
	psdx = (1/(Fs*N)) * np.square(abs(xdft))
	psdx[1:len(psdx)-1] = 2*psdx[1:len(psdx)-1]
	#print(Fs)
	#print(Fs/N)
	#print(N)
	freq = np.arange(0, Fs/2, Fs/N)
	return psdx, freq




# apply the PIO model to the data WMLES dataset
ii=1
jj=1
cali_size = tau_w_star.shape
#print(cali_size)
#cali_size = [4348, 10, 2]
tau_wp = np.zeros(cali_size, dtype=np.float64)
folder_bin = folder + 'alpha.bin'
f = open(folder_bin, 'rb')
alpha = np.fromfile(f, dtype=np.float64)
psdx_size = [cali_size[0]/2, cali_size[1], cali_size[2]]
tau_wp_psdx = np.zeros(psdx_size, dtype=np.float64)
tau_wp_f = np.zeros(psdx_size, dtype=np.float64)
print(alpha[0], alpha[1])

for j in range(0, cali_size[1]):
	for k in range(0, cali_size[2]):
		tau_wp[:,j,k], Fs = PIO_model_apply(u_WMLES[:,ii,jj], tau_w_star[:,j,k],t_WMLES, t_tau[:, k], d, retau, alpha[k])
		tau_wp_psdx[:,j,k], tau_wp_f[:,j,k] = psdx(tau_wp[:,j,k], Fs)	 
print(rms(tau_wp[:,0,0]))
#print(tau_wp.shape)
print('Model applied and psdx computed')

#Mean PSDX
print(tau_wp_psdx.shape)
summ = np.zeros(psdx_size[0], dtype=np.float64)
for j in range(0, psdx_size[1]):
	for k in range(0, psdx_size[2]):
		summ = summ + tau_wp_psdx[:,j,k]
tau_wp_psdx_mean = summ/(float(j*k))

# Plots
print('Begin plotting')
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
# PSDX spectra of shear stress flucctuations
fig1 = plt.figure()
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
ax = fig1.add_subplot(111)
psdx_plot = np.zeros(tau_wp_psdx_mean.shape, dtype=np.float64)
for i in range(0, len(tau_wp_psdx_mean)):
	psdx_plot[i] = (tau_wp_f[i,0,0] * tau_wp_psdx_mean[i])/1.035**2
plt.plot(2*np.pi*tau_wp_f[:,0,0], psdx_plot, 'b-')
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylim([10**-4, 10**0])
plt.xlim([10**-3, 2**1])
#plt.xlabel(r"$w^{+}$", fontsize=16)
#plt.ylabel(r"$w^{+} E(w^{+})/ \tau_{w}^{+2}$",fontsize=16)
plt.savefig('test_fig_2.eps', format='eps', dpi=1000)
print('Plots saved')



# Output
sys.stderr.close()
sys.stderr = sys.__stderr__




