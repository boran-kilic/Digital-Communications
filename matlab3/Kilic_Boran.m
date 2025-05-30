clear; close all; clc;
%% Question 1
%% part a
rng(37);

Fs = 5000;
T = 0.1;
t = 0:1/Fs:T-1/Fs;
N = length(t);
h = ones(1,N);

s1 =cos(2*pi*250*t).*h;
s2 =cos(2*pi*500*t).*h;

numSymbols = 5;
tx_symbols = randi(4, 1, numSymbols);            
x = zeros(1,numSymbols*N);

for i = 1:numSymbols
    sym = tx_symbols(i); 
    if sym == 1         %"00"
        nextsym = s1;   
    elseif sym == 2     %"01"
        nextsym = -s1;  
    elseif sym == 3     %"10"
        nextsym = s2;
    elseif sym == 4     %"11"
        nextsym = -s2;
    end
    x((i-1)*N+1:(i)*N) = nextsym;
end

time = 0:1/Fs:T*numSymbols-1/Fs;
plot(time, x);
xlabel('Time (s)');
ylabel('x(t)');
title('Modulated Signal x(t)');
grid on;

%% part b

phi1 = 2*sqrt(5) *cos(2*pi*250*t);
phi2 = 2*sqrt(5) *cos(2*pi*500*t);

figure;
subplot(2,1,1);
plot(t, phi1, 'b');
title('$\psi_1(t) = 2\sqrt{5} \cos(250\pi t)$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, phi2, 'r');
title('$\psi_2(t) = 2\sqrt{5} \cos(500\pi t)$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% part c

noise_variances = [1e-4, 1e-2, 1e0];

for i = 1:length(noise_variances)
    noise = sqrt(noise_variances(i)) * randn(size(x));
    x_noisy = x + noise; 
    
    figure; 
    plot(time, x, 'LineWidth', 1.2); hold on; 
    plot(time, x_noisy); 
    title(['Original and Noisy Signal (\sigma^2 = ' num2str(noise_variances(i)) ')']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original', 'Noisy');
    grid on;
end

%% part e

numSymbols = 1e5;
tx_symbols = randi(4, 1, numSymbols);            
x = zeros(1,numSymbols*N);

for i = 1:numSymbols
    sym = tx_symbols(i); 
    if sym == 1         %"00"
        nextsym = s1;   
    elseif sym == 2     %"01"
        nextsym = -s1;  
    elseif sym == 3     %"10"
        nextsym = s2;
    elseif sym == 4     %"11"
        nextsym = -s2;
    end
    x((i-1)*N+1:(i)*N) = nextsym;
end

SNR_dB = 0:1:12;
SNR = 10.^(SNR_dB/10);

Es =  trapz(t,s1.^2);  
N0 = Es ./ SNR;
qfunc = @(x) 0.5 * erfc(x / sqrt(2));
Pe_theory = 2*qfunc(sqrt(Es ./ N0)) - qfunc(sqrt(Es ./ N0)).^2;

Pe_sim = zeros(size(SNR));

for k = 1:length(SNR)
    noise = sqrt(N0(k)/2*Fs) * randn(size(x));
    r = x + noise;

    rx_symbols = zeros(1, numSymbols);
    for i = 1:numSymbols
        segment = r((i-1)*N + (1:N));
        r1 = sum(segment .* phi1);
        r2 = sum(segment .* phi2);

        if r1 >= r2 && r1 >= -r2
            rx_symbols(i) = 1;      %"00"
        elseif r1 < r2 && r1 <= -r2
            rx_symbols(i) = 2;      %"01"
        elseif r1 < r2 && r1 > -r2
            rx_symbols(i) = 3;      %"10"
        elseif r1 >= r2 && r1 < -r2
            rx_symbols(i) = 4;      %"11"
        end

    end

    Pe_sim(k) = sum(rx_symbols ~= tx_symbols) / numSymbols;
end

figure;
semilogy(SNR_dB, Pe_sim, 'bo-', 'LineWidth', 1.5);
hold on;
semilogy(SNR_dB, Pe_theory, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Symbol Error Probability');
legend('Simulated', 'Theoretical');
title('Symbol Error Probability vs. SNR for Coherent FSK');

%% Question 2
%% part a
clear;
rng(37);

Fs = 1000;
T = 0.01;
t = 0:1/Fs:T-1/Fs;
N = length(t);
h = ones(1,N);

s =sin(2*pi*100*t).*h;

numSymbols = 5;

tx_symbols = randi(2, 1, numSymbols);            
x = zeros(1,numSymbols*N);

for i = 1:numSymbols
    sym = tx_symbols(i); 
    if sym == 1         %"0"
         nextsym = s;
    elseif sym == 2     %"1"
         nextsym = -s;
    end
    x((i-1)*N+1:(i)*N) = nextsym;
end

time = 0:1/Fs:T*numSymbols-1/Fs;
figure;
plot(time, x);
xlabel('Time (s)');
ylabel('x(t)');
title('Modulated Signal x(t)');
grid on;

%% part b

phi1 = 10*sqrt(2) *sin(2*pi*100*t);

figure;
plot(t, phi1, 'b');
title('$\psi_1(t) = 10\sqrt{2} \sin(100\pi t)$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% part c

noise_variances = [1e-4, 1e-2, 1e0];

for i = 1:length(noise_variances)
    noise = sqrt(noise_variances(i)) * randn(size(x));
    x_noisy = x + noise; 
    
    figure; 
    plot(time, x, 'LineWidth', 1.2); hold on; 
    plot(time, x_noisy); 
    title(['Original and Noisy Signal (\sigma^2 = ' num2str(noise_variances(i)) ')']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original', 'Noisy');
    grid on;
end

%% part e 

numSymbols = 1e5;
tx_symbols = randi(2, 1, numSymbols)-1;            
x = zeros(1,numSymbols*N);

for i = 1:numSymbols
    sym = tx_symbols(i); 
    if sym == 0         %"0"
        nextsym = s;   
    elseif sym == 1     %"1"
        nextsym = -s;  
    end
    x((i-1)*N+1:(i)*N) = nextsym;
end

SNR_dB = 0:1:9;
SNR = 10.^(SNR_dB/10);

Es =  trapz(t,s.^2);  
N0 = Es ./ SNR;
qfunc = @(x) 0.5 * erfc(x / sqrt(2));
Pe_theory = qfunc(sqrt(2*Es ./ N0));

Pe_sim = zeros(size(SNR));

for k = 1:length(SNR)
    noise = sqrt(N0(k)/2*Fs) * randn(size(x));
    r = x + noise;

    rx_symbols = zeros(1, numSymbols);
    for i = 1:numSymbols
        segment = r((i-1)*N + (1:N));
        r1 = sum(segment .* phi1);

        if r1 > 0
            rx_symbols(i) = 0;      %"0"
        elseif r1 < 0
            rx_symbols(i) = 1;      %"1"
        end

    end

    Pe_sim(k) = sum(rx_symbols ~= tx_symbols) / numSymbols;
end

figure;
semilogy(SNR_dB, Pe_sim, 'bo-', 'LineWidth', 1.5);
hold on;
semilogy(SNR_dB, Pe_theory, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Symbol Error Probability');
legend('Simulated', 'Theoretical');
title('Symbol Error Probability vs. SNR for BPAM with equal priors');

%% part g

N0 = 2* 1e-2;

alphas = linspace(1e-4, 0.5 - 1e-4, 20);

Pe_sim1 = zeros(size(alphas));
Pe_sim2 = zeros(size(alphas));

for k = 1:length(alphas)
    tx_symbols = zeros(1, numSymbols); 
    tx_symbols(rand(1, numSymbols) < alphas(k)) = 1;

    x = zeros(1,numSymbols*N);
    
    for i = 1:numSymbols
        sym = tx_symbols(i); 
        if sym == 0         %"0"
             nextsym = s;
        elseif sym == 1     %"1"
             nextsym = -s;
        end
        x((i-1)*N+1:(i)*N) = nextsym;
    end
    
    noise = sqrt((N0/2)*Fs) * randn(size(x));
    r = x + noise;
    r_th = log(alphas(k)/(1-alphas(k)))*(N0*Fs/(4*sqrt(Es)));
    
    rx_symbols = zeros(1, numSymbols);
    for i = 1:numSymbols
        segment = r((i-1)*N + (1:N));
        r1 = sum(segment .* phi1);

        if r1 > r_th
            rx_symbols(i) = 0;      %"0"
        elseif r1 < r_th
            rx_symbols(i) = 1;      %"1"
        end

    end

    Pe_sim1(k) = sum(rx_symbols ~= tx_symbols) / numSymbols;

    rx_symbols = zeros(1, numSymbols);
    for i = 1:numSymbols
        segment = r((i-1)*N + (1:N));
        r1 = sum(segment .* phi1);

        if r1 > 0
            rx_symbols(i) = 0;      %"0"
        elseif r1 < 0
            rx_symbols(i) = 1;      %"1"
        end

    end

    Pe_sim2(k) = sum(rx_symbols ~= tx_symbols) / numSymbols;
end


figure;
stem(alphas, Pe_sim1, 'bo-', 'LineWidth', 1.5); hold on;
stem(alphas, Pe_sim2, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('\alpha');
ylabel('Symbol Error Probability');
legend('MAP receiver','MLE receiver');
title('Symbol Error Probability vs \alpha');
