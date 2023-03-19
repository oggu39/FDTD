import numpy as np
import matplotlib.pyplot as plt


E_reflect = np.loadtxt('exported_data\\reflected_field.txt');
t = np.loadtxt('exported_data\\t_data.txt');
E_transmit = np.loadtxt('exported_data\\transmitted_field.txt');
E_src = np.loadtxt('exported_data\\source_field.txt');


E_reflect[t<(2.85*(10**(-9)))] = 0; # Remove incident initial wave since its not reflected

plt.figure(1);
plt.plot(t,E_reflect, label = "Reflected field");
plt.plot(t,E_transmit, label = "Transmitted field");
plt.plot(t,E_src, label = "Source field");
plt.legend();
plt.xlabel('Simulation time[s]');
plt.ylabel('Field strength[V/m]');
plt.title("Time domain");

dt = t[1]-t[0];
fs = 1/dt;

plt.figure(2);
F_reflect = np.abs(np.fft.fft(E_reflect))
F_transmit = np.abs(np.fft.fft(E_transmit)); # Fourier transform
F_src = np.abs(np.fft.fft(E_src));
freq = fs*np.fft.fftfreq(t.shape[-1]);
f_max = 5*10**9;

F_reflect = F_reflect[(freq>(-f_max))*(freq<f_max)];
F_transmit = F_transmit[(freq>(-f_max))*(freq<f_max)];
F_src = F_src[(freq>(-f_max))*(freq<f_max)];
freq = freq[(freq>(-f_max))*(freq<f_max)];


norm_factor = max(F_src);
F_src = F_src/norm_factor;
F_reflect = ((F_reflect/norm_factor)/F_src)**2;
F_transmit = ((F_transmit/norm_factor)/F_src)**2;
total = F_reflect+F_transmit;
print("Mean reflection coefficient for glass: " + str(np.mean(F_reflect[0:-2])));


plt.plot(freq[0:-2], (F_reflect[0:-2]), label = 'Reflection');
plt.plot(freq[0:-2], (F_transmit[0:-2]),label = 'Transmission');
plt.plot(freq[0:-2], total[0:-2], label = 'Total');




plt.legend();
plt.xlabel('Frequency');
plt.ylabel('Transform field');
plt.title('Frequency domain');

plt.show();
