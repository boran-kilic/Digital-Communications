clear
close all
clc
%%
%part a

z = linspace(-10, 10, 1000);
fZ = (1/2) * exp(-abs(z));
figure
plot(z, fZ);
xlabel('z'); ylabel('f_Z(z)');
title('PDF of Z = X - Y');
grid on;

%%
%part b

rng(17);
N1 = 100;    
N2 = 100000; 

X1 = exprnd(1, N1, 1);  
Y1 = exprnd(1, N1, 1);  
Z1 = X1 - Y1;

X2 = exprnd(1, N2, 1); 
Y2 = exprnd(1, N2, 1);  
Z2 = X2 - Y2;

figure
subplot(1, 2, 1);
histogram(Z1, 'Normalization', 'pdf', 'EdgeColor', 'none', 'BinWidth', 0.1);
hold on;
plot(z, fZ, 'r', 'LineWidth', 2);
title('Histogram of Z for 100 realizations');
xlabel('Z');
ylabel('Density');
legend('Empirical', 'Theoretical');
grid on;

subplot(1, 2, 2);
histogram(Z2, 'Normalization', 'pdf', 'EdgeColor', 'none', 'BinWidth', 0.1);
hold on;
plot(z, fZ, 'r', 'LineWidth', 2);
title('Histogram of Z for 100,000 realizations');
xlabel('Z');
ylabel('Density');
legend('Empirical', 'Theoretical');
grid on;

%%
%part c
% Z = Z2;
Z = Z2(Z2 >= -10 & Z2 <= 10);
N = 100000;

empirical_avg_Z2 = mean(Z.^2);
fprintf('Empirical average of Z^2: %.4f\n', empirical_avg_Z2);
fprintf('\n');

%%
%part d

num_levels = 8;
Z_min = min(Z);
Z_max = max(Z);
delta = (Z_max - Z_min) / num_levels ;

quantized_Z = Z;

for i = 1:length(Z)
    if Z(i) == Z_max
    quantized_Z(i) = Z_min + (delta * floor((Z(i) - Z_min) / delta)) - delta/2;
    else
    quantized_Z(i) = Z_min + (delta * floor((Z(i) - Z_min) / delta)) + delta/2;
    end
end

figure
plot(Z, quantized_Z, 'r.', 'LineWidth', 1);
xlabel('$z$', 'Interpreter', 'latex');
ylabel('$Q(z)$', 'Interpreter', 'latex');
title('Quantization Regions and Reconstruction Levels');
grid on;

quantization_error = Z - quantized_Z;

P_e = mean(quantization_error.^2);

P_signal = mean(Z.^2);

SQNR_dB = 10 * log10(P_signal / P_e);

fprintf('Uniform quantizer\n');
fprintf('Average power of the quantization error: %.4f\n', P_e);
fprintf('Signal power: %.4f\n', P_signal);
fprintf('SQNR (dB): %.4f\n', SQNR_dB);
fprintf('\n');

%%
%part e

quantization_levels = Z_min + delta/2 : delta : Z_max - delta/2; 
EDGES = Z_min:delta: Z_max;

fprintf('Quantization Boundaries: \n');
disp(EDGES);
fprintf('Reconstruction Levels: \n');
disp(quantization_levels);

freq = histcounts(quantized_Z, EDGES); 

pmf = freq / N;

figure;
bar(quantization_levels, pmf, 'FaceColor', 'b');
xlabel('$\tilde{Z}$', 'Interpreter', 'latex');
ylabel('Probability');
title('Probability Mass Function (PMF) of Quantized Z');
grid on;

%%
%part f

F_Z = @(z) (z < 0) .* (1/2 * exp(-abs(z))) + (z >= 0) .* (1 - 1/2 * exp(-abs(z)));
F_Z_inv = @(p) ((0 < p) & (p <= 0.5)) .* log(2*p) + ((1 > p) & (p > 0.5)) .* (-log(2*(1 - p)));


z_vals = linspace(-5, 5, 1000);
F_vals = F_Z(z_vals);
figure;
plot(z_vals, F_vals, 'LineWidth', 2);
xlabel('$z$', 'Interpreter', 'latex');
ylabel('$F_Z(z)$', 'Interpreter', 'latex');
title('CDF of Z');
grid on;

p_vals = linspace(0.01, 0.99, 1000);  
F_Z_inv_vals = F_Z_inv(p_vals);
figure;
plot(p_vals, F_Z_inv_vals, 'LineWidth', 2);
xlabel('$p$', 'Interpreter', 'latex');
ylabel('$F_Z^{-1}(p)$', 'Interpreter', 'latex');
title('Inverse CDF of Z');
grid on;

%%
%part g

compressed_Z = F_Z(Z);

num_levels = 8;
min_compressed = min(compressed_Z);
max_compressed = max(compressed_Z);
delta = (max_compressed - min_compressed) / num_levels;

quantized_compressed_Z = compressed_Z;
for i = 1:length(compressed_Z)
    if compressed_Z(i) == max_compressed
        quantized_compressed_Z(i) = min_compressed + (delta * floor((compressed_Z(i) - min_compressed) / delta)) - delta / 2;
    else
        quantized_compressed_Z(i) = min_compressed + (delta * floor((compressed_Z(i) - min_compressed) / delta)) + delta / 2;
    end
end

reconstructed_Z = F_Z_inv(quantized_compressed_Z);

figure;
plot(Z, reconstructed_Z, 'r.', 'LineWidth', 1);
xlabel('$z$', 'Interpreter', 'latex');
ylabel('$Q(z)$', 'Interpreter', 'latex');
title('Quantization Regions and Reconstruction Levels');
grid on;


%%
%part h
quantization_error2 = Z - reconstructed_Z;

P_e = mean(quantization_error2.^2);

P_signal = mean(Z.^2);

SQNR_dB = 10 * log10(P_signal / P_e);

fprintf('Non-uniform PCM quantizer\n');
fprintf('Average power of the quantization error: %.4f\n', P_e);
fprintf('Signal power: %.4f\n', P_signal);
fprintf('SQNR (dB): %.4f\n', SQNR_dB);
fprintf('\n');

%%
%part i

quantization_levels = min_compressed + delta/2 : delta : max_compressed - delta/2; 
EDGES = min_compressed:delta:max_compressed;

fprintf('Quantization Boundaries: \n');
disp(EDGES);
fprintf('Reconstruction Levels: \n');
disp(quantization_levels);

freq = histcounts(quantized_compressed_Z, EDGES);

pmf = freq/N;

figure;
bar(quantization_levels, pmf, 'FaceColor', 'b');
xlabel('$\tilde{Z}$', 'Interpreter', 'latex');  % Quantized values
ylabel('Probability', 'Interpreter', 'latex');
title('Probability Mass Function (PMF) of Quantized Z');
grid on;

%%
%part j

num_levels = 8;

Z_min = min(Z);
Z_max = max(Z);

%Step 1 Initilize
boundary = linspace(Z_min, Z_max, num_levels + 1); 
reconstruction_levels = (boundary(1:end-1) + boundary(2:end)) / 2;

max_iter = 100;
epsilon = 1e-6;

for iter = 1:max_iter
    % Step 2 E-step
    quantized_Z = zeros(length(Z), 1);
    for i = 1:length(Z)
        [~, idx] = min(abs(Z(i) - reconstruction_levels));
        quantized_Z(i) = reconstruction_levels(idx);
    end
    
    
    new_reconstruction_levels = zeros(1, num_levels);
    for k = 1:num_levels
        region_values = Z(quantized_Z == reconstruction_levels(k));
        if ~isempty(region_values)
            new_reconstruction_levels(k) = mean(region_values);  
        else
            new_reconstruction_levels(k) = reconstruction_levels(k);  
        end
    end

    % Step 3 M-step
    new_boundary = boundary;
    for k = 1:num_levels-1
        new_boundary(k+1) = (new_reconstruction_levels(k) + new_reconstruction_levels(k+1)) / 2;
    end
    
    % Step 4 Convergence
    if max(abs(new_reconstruction_levels - reconstruction_levels)) < epsilon
        break;
    end
    
    reconstruction_levels = new_reconstruction_levels;
    boundary = new_boundary;
end


figure;
hold on;
for k = 1:num_levels
    plot([boundary(k), boundary(k+1)], [reconstruction_levels(k), reconstruction_levels(k)], 'r-', 'LineWidth', 3);  % Quantization regions
end
xlabel('$z$', 'Interpreter', 'latex');
ylabel('$Q(z)$', 'Interpreter', 'latex');
title('Quantization Regions and Reconstruction Levels');
grid on;

%%
%part k

quantized_Z = zeros(length(Z), 1);
for i = 1:length(Z)
    [~, idx] = min(abs(Z(i) - reconstruction_levels));
    quantized_Z(i) = reconstruction_levels(idx);
end

quantization_error = Z - quantized_Z;

P_e = mean(quantization_error.^2);

P_signal = mean(Z.^2);

SQNR_dB = 10 * log10(P_signal / P_e);

fprintf('LLyod-Max Quantizier\n');
fprintf('Average power of the quantization error: %.4f\n', P_e);
fprintf('Signal power: %.4f\n', P_signal);
fprintf('SQNR (dB): %.4f\n', SQNR_dB);
fprintf('\n');

fprintf('Quantization Boundaries: \n');
disp(boundary);
fprintf('Reconstruction Levels: \n');
disp(reconstruction_levels);

%%
%part l

freq = histcounts(quantized_Z, boundary);
pmf = freq/N;
figure;
bar(reconstruction_levels, pmf, 'FaceColor', 'b');
xlabel('$\tilde{Z}$', 'Interpreter', 'latex');  % Quantized values
ylabel('Probability', 'Interpreter', 'latex');
title('Probability Mass Function (PMF) of Quantized Z');
grid on;
