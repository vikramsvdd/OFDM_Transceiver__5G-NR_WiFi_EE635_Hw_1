% Clear workspace
clear all;

% System parameters
num_subcarriers = 8;           % Number of subcarriers
cyclic_prefix_length = 4;      % Length of cyclic prefix
noise_variance = 0.2;          % Noise variance
num_simulations = 500;        % Number of Monte Carlo simulations

% Channel impulse response
channel_impulse_response = [0.5 + 0.2j, -0.6 + 0.1j, 0.2 - 0.25j].';
channel_frequency_response = fft(channel_impulse_response, num_subcarriers);

% Input data and modulation
input_bits = [1 0 1 1 0 0 0 1];  % Input bit sequence
bpsk_symbols = 2 * input_bits - 1;  % BPSK modulation: 0 -> -1, 1 -> 1

% Hadamard transform and IFFT for waveform generation
hadamard_matrix = hadamard(num_subcarriers);
frequency_domain_signal = hadamard_matrix * bpsk_symbols.';
time_domain_signal = ifft(frequency_domain_signal);

% Add cyclic prefix
transmitted_signal = [time_domain_signal(end-cyclic_prefix_length+1:end); time_domain_signal];

% Initialize error counter
total_bit_errors = 0;

% Monte Carlo simulations
for sim_index = 1:num_simulations
    % Channel convolution and noise addition
    noise = sqrt(noise_variance/2) * (randn(length(transmitted_signal) + length(channel_impulse_response) - 1, 1) + ...
                                      1j * randn(length(transmitted_signal) + length(channel_impulse_response) - 1, 1));
    received_signal = conv(channel_impulse_response, transmitted_signal) + noise;
    
    % Remove cyclic prefix and extra samples
    useful_signal = received_signal(cyclic_prefix_length + 1 : num_subcarriers + cyclic_prefix_length);
    
    % FFT and Zero-Forcing Equalization
    frequency_domain_received = fft(useful_signal);
    equalized_signal = frequency_domain_received ./ channel_frequency_response;
    
    % Inverse Hadamard transform
    estimated_symbols = (1/num_subcarriers) * hadamard_matrix * equalized_signal;
    
    % BPSK demodulation using minimum distance detection
    demodulated_bits = real(estimated_symbols) > 0;
    
    % Count bit errors
    bit_errors = sum(demodulated_bits ~= input_bits.');
    total_bit_errors = total_bit_errors + bit_errors;
end

% Calculate average Bit Error Rate (BER)
average_ber = total_bit_errors / (num_simulations * num_subcarriers);

% Display results
fprintf('Average Bit Error Rate: %e\n', average_ber);
