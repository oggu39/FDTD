clear
close all
clc

reflection_file = "reflected_field.txt";
transmission_file = "transmitted_field.txt";
source_file = "source_field";
t_file = "t_data.txt";

figure()
subplot(1,2,1);
[E_r, E_t,E_src, t] = read_data(reflection_file, transmission_file, source_file,  t_file);
E_r(t<0.25e-8) = 0;
dt = t(2)-t(1);
fs = 1/dt;
f = fs*linspace(-1/2,1/2,length(t));
plot(t,E_r, t, E_t, t, E_src);
xlabel('Time[s]');
ylabel('Field strength[V/m]');
legend('Reflected field', 'Transmitted field', 'Source field');
axis([0 1.2e-8 -2 2])


F_r = fftshift(fft(E_r));
F_t = fftshift(fft(E_t));
F_src = fftshift(fft(E_src));
subplot(1,2,2);
plot(f,abs(F_src),f,abs(F_r),f, abs(F_t));
xlabel('Frequency');
ylabel('Fourier transform');
legend('Source', 'Reflected', 'Transmitted');
a = 0.05;
axis([-a*max(f) a*max(f) 0 max(abs(F_src))]);

