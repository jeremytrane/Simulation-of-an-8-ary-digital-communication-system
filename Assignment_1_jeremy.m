clear;
clc;

% Parameters
d = 1; 
constellation = d/2 * [1 + 3i, -1 + 3i, 1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i, 1 - 3i, -1 - 3i];
normalized_constellation = constellation / sqrt(1.5);
Eb_N0_dB = -6:2:12;
Es_N0_dB = Eb_N0_dB + 10*log10(log2(8));
no_symbols_sent = 100000;

% Generate data and corresponding signals
data = randi([0, 7], no_symbols_sent, 1);
signal_without_noise = normalized_constellation(data + 1);

SER_simulated = zeros(1, length(Eb_N0_dB));
SER_theoretical = zeros(1, length(Eb_N0_dB));
BER = zeros(1, length(Eb_N0_dB));

for idx = 1:length(Eb_N0_dB)
    noiseSigma = sqrt(1 / (2 * log2(8) * 10^(Eb_N0_dB(idx) / 10)));
    noise = noiseSigma * (randn(size(signal_without_noise)) + 1i*randn(size(signal_without_noise)));
    received = signal_without_noise + noise;
    
    % Distance calculation using meshgrid and reshaping
    [R, C] = meshgrid(received, normalized_constellation);
    distance = abs(R - C).^2;
    [~, minDist] = min(distance, [], 1);
    
    SER_simulated(idx) = sum(minDist' - 1 ~= data) / no_symbols_sent;
    SER_theoretical(idx) = (5/2)*qfunc(sqrt(10^(Es_N0_dB(idx)/10)/3)) - (3/2)*qfunc(sqrt(10^(Es_N0_dB(idx)/10)/3))^2;
    BER(idx) = 1/log2(8) * SER_simulated(idx);
end

% Plotting
figure();
semilogy(Es_N0_dB, SER_theoretical, 'r-');
hold on;
semilogy(Es_N0_dB, SER_simulated, 'b-');
legend('Theory', 'Simulated');
xlabel('Es/N0(dB)');
ylabel('SER');
title('Comparison of theoretical vs simulated SER');
grid on;

figure();
semilogy(Eb_N0_dB, BER, 'r-*');
hold on;
xlabel('Eb/N0(dB)');
ylabel('BER');
title('Plot of BER');
grid on;
