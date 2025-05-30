clear; close all; clc;
%% Question 1
%% part b
rng(29)
N = 1e6;                                 
t = 10000:11000;
W = zeros(1,N); % Random Process W(t0, n), actually a rv   
V = randn(1,N); % Gaussian noise 

%%
for n = 2:N
    W(n) = 0.9*W(n-1) + V(n);
end

%%
figure;
plot(t,W(t));
title("Realizations of Wn");
xlabel("n");ylabel("Wn");

%% part c

Y = zeros(1,N);
for n = 1:N
    for m = 1:min(n,20)   
        Y(n) = Y(n) + (1/20)*W(n-m+1);
    end
end

%%
figure;
plot(t,W(t));
title("Realizations of Wn and Yn");
hold on
plot(t,Y(t));
legend ("W_n","Y_n")
xlabel("n");ylabel("Value");

%% part d

N = 1e6;
window = hamming(1024);
noverlap = 512;
nfft = 2048;              

[pxxWn, fWn] = pwelch(W, window, noverlap, nfft);
[pxxYn, fYn] = pwelch(Y, window, noverlap, nfft);

H_estimated = sqrt(pxxYn ./ pxxWn);  

n = 0:19;
h = ones(1,20)/20;                 
[H_analytical, fH] = freqz(h, 1, nfft, 2*max(fWn));  

figure;
plot(fWn, 10*log10(pxxWn));
xlim([min(fWn) max(fWn)]);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
hold on
plot(fYn, 10*log10(pxxYn));
xlim([min(fYn) max(fYn)]);
title('PSD estimates of Wn and Yn');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend ("S_W","S_Y")
grid on;
hold off

figure;
plot(fWn, H_estimated, 'b', 'DisplayName', 'Estimated');
hold on;
plot(fH, abs(H_analytical), 'r--', 'DisplayName', 'Analytical');
xlim([min(fH) max(fH)]);
title('Magnitude Response Comparison');
xlabel('Frequency (Hz)');
ylabel('|H(f)|');
legend;
grid on;
hold off

%% part e

Y = zeros(1,N);
for n = 1:N
    for m = 1:min(n,20)   
        Y(n) = Y(n) + (1/20)*V(n-m+1);
    end
end

%%
figure;
plot(t,V(t));
title("Realizations of Vn and Yn");
hold on
plot(t,Y(t));
legend ("Vn","Yn")
xlabel("n");ylabel("Value");

%%
N = 1e6;
window = hamming(1024);
noverlap = 512;
nfft = 2048;              

[pxxVn, fVn] = pwelch(V, window, noverlap, nfft);
[pxxYn, fYn] = pwelch(Y, window, noverlap, nfft);

H_estimated = sqrt(pxxYn ./ pxxVn);  

n = 0:19;
h = ones(1,20)/20;                 
[H_analytical, fH] = freqz(h, 1, nfft, 2*max(fVn));  

figure;
plot(fVn, 10*log10(pxxVn));
xlim([min(fVn) max(fVn)]);

hold on
plot(fYn, 10*log10(pxxYn));
xlim([min(fYn) max(fYn)]);

xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD estimates of V_n and Y_n');
legend ("S_V","S_Y")
grid on;
hold off

figure;
plot(fWn, 10*log10(pxxWn));
xlim([min(fWn) max(fWn)]);

hold on
plot(fVn, 10*log10(pxxVn));
xlim([min(fVn) max(fVn)]);

xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD estimates of V_n and W_n');
legend ("S_W","S_V")
grid on;
hold off

figure;
plot(fVn, H_estimated, 'b', 'DisplayName', 'Estimated');
hold on;
plot(fH, abs(H_analytical), 'r--', 'DisplayName', 'Analytical');
axis tight
title('Magnitude Response Comparison');
xlabel('Frequency (Hz)');
ylabel('|H(f)|');
legend;
grid on;
hold off


%% Question 2
%% part a
clear;

T = 0.01;
Fs = 100000;
t = [0:T*Fs]'/Fs;
mt = sin(200*pi*t);
figure
plot(t',mt);
title('Message Signal');
xlabel('Time (s)');
ylabel('Amplitude');

fc = 5000;
Ac1=1;
xt = Ac1.*mt.*cos(2*pi.*(fc).*t);
figure
plot(t,xt);
title('Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

%% part b
Ac2 = 2;
delta = [0 50 100 500 1000];

[b, a] = butter(6, 2*2500/Fs);

figure;
hold on;
colors = lines(length(delta));

for k = 1:length(delta)
    x_local = Ac2 * cos(2*pi*(fc + delta(k)) * t);
    vt = xt .* x_local;
    v0t = filtfilt(b, a, vt);
    plot(t, v0t, 'DisplayName', ['\Delta = ' num2str(delta(k))], ...
         'Color', colors(k,:));
end

title('Demodulated Signals for Different \Delta');
xlabel('Time (s)');
ylabel('Amplitude');
legend show;
hold off;

%% part d
Ac2 = 2;
phi = [0 pi/6 pi/3 pi/2];

figure;
hold on;
colors = lines(length(phi));

for k = 1:length(phi)
    x_local = Ac2 * cos(2*pi*fc.*t + phi(k));
    vt = xt .* x_local;
    v0t = filtfilt(b, a, vt);
    plot(t, v0t, 'DisplayName', ['\phi = ' num2str(delta(k))], ...
         'Color', colors(k,:));
end

title('Demodulated Signals for Different \phi');
xlabel('Time (s)');
ylabel('Amplitude');
legend show;
hold off;
