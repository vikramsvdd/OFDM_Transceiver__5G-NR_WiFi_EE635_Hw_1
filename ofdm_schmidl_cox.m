% Simulate an optimized Schmidl-Cox synchronization algorithm for an OFDM system
% with 64 subcarriers, 10 MHz bandwidth, 1.74 microseconds timing offset,
% and an AWGN channel with 6 dB SNR.

%This system is actively used in 5G-NR/WiFi and was used for LTE

% System parameters
SUBCARRIERS = 64;          % Number of OFDM subcarriers
BANDWIDTH = 1e7;           % Bandwidth in Hz
TIMING_OFFSET = 1.74e-6;   % Timing offset in seconds
CYCLIC_PREFIX = 5;         % Length of cyclic prefix
SUBCARRIER_SPACING = BANDWIDTH / SUBCARRIERS;  % Subcarrier spacing (Hz)
SYMBOL_DURATION = 1 / SUBCARRIER_SPACING;      % Symbol duration (seconds)
FREQ_OFFSET = 200e3 / SUBCARRIER_SPACING;      % Normalized frequency offset

% Generate PN sequences for integer frequency estimation
%Here we are generating two sequences, one is for odd subcarriers, the
%other one is for even
pn_gen1 = comm.PNSequence('Polynomial', [4 1 0], 'SamplesPerFrame', SUBCARRIERS/2, 'InitialConditions', [0 1 0 0]);
pn_gen2 = comm.PNSequence('Polynomial', [4 3 0], 'SamplesPerFrame', SUBCARRIERS/2, 'InitialConditions', [0 0 1 0]);

% Generate OFDM signal with cyclic prefix
TOTAL_SYMBOLS = 300;
tx_symbols = zeros(SUBCARRIERS, TOTAL_SYMBOLS);
tx_waveform = zeros(SUBCARRIERS, TOTAL_SYMBOLS);
tx_waveform_cp = zeros(SUBCARRIERS + CYCLIC_PREFIX, TOTAL_SYMBOLS);

for sym_idx = 1:TOTAL_SYMBOLS  %OFDM symbol index, 3 categories
    if sym_idx ~= 2
        % Regular data symbols
        tx_symbols(:, sym_idx) = (2 * round(rand(SUBCARRIERS, 1)) - 1) * (0.707 + 1i*0.707);
    else
        % Frequency synchronization symbol
        tx_symbols(1:2:end, sym_idx) = (2 * pn_gen1() - 1) * (0.707 + 1i*0.707);
        tx_symbols(2:2:end, sym_idx) = (2 * pn_gen2() - 1) * (0.707 + 1i*0.707);
    end
    
    % Apply IFFT and frequency offset, using ifft not traditional frequency
    % generator based OFDM
    tx_waveform(:, sym_idx) = ifft(tx_symbols(:, sym_idx)) .* exp(1i * 2 * pi * FREQ_OFFSET * (0:SUBCARRIERS-1)' / SUBCARRIERS);
    
    % Special processing for the first symbol (synchronization symbol)
    if sym_idx == 1
        tx_waveform(SUBCARRIERS/2 + 1:end, sym_idx) = tx_waveform(1:SUBCARRIERS/2, sym_idx) * exp(1i * pi * FREQ_OFFSET);
    end
    
    % Add cyclic prefix
    tx_waveform_cp(:, sym_idx) = [tx_waveform(end-CYCLIC_PREFIX+1:end, sym_idx); tx_waveform(:, sym_idx)];
end

% Serialize the transmit waveform
tx_waveform_serial = reshape(tx_waveform_cp, 1, []);

% Simulation parameters
NUM_ITERATIONS = 500;
SNR_DB = 6;
timing_errors = zeros(1, NUM_ITERATIONS);

% Main simulation loop
for iter = 1:NUM_ITERATIONS
    % Add AWGN noise
    signal_power = mean(abs(tx_waveform_serial).^2);
    noise_power = signal_power / 10^(SNR_DB/10);
    noise = sqrt(noise_power/2) * (randn(size(tx_waveform_serial)) + 1i*randn(size(tx_waveform_serial)));
    rx_waveform = tx_waveform_serial + noise;
    
    % Add delay
    sample_delay = round(TIMING_OFFSET * BANDWIDTH);
    rx_waveform_delayed = [zeros(1, sample_delay), rx_waveform(1:end-sample_delay)];
    
    % Schmidl-Cox autocorrelation - the core part lives here
    correlation_window = 200;
    autocorr_func = zeros(1, correlation_window);
    for window_idx = 1:correlation_window
        rx_segment = rx_waveform_delayed(window_idx + CYCLIC_PREFIX : window_idx + SUBCARRIERS + CYCLIC_PREFIX - 1);
        autocorr_func(window_idx) = sum(rx_segment(SUBCARRIERS/2+1:end) .* conj(rx_segment(1:SUBCARRIERS/2))) / ...
                                    sum(abs(rx_segment(SUBCARRIERS/2+1:end)).^2);
    end
    
    % Find symbol timing
    [~, estimated_start] = max(abs(autocorr_func));
    true_delay = sample_delay + 1;
    timing_errors(iter) = (true_delay - estimated_start) / BANDWIDTH;
    
    % Plot autocorrelation function for the first iteration
    if iter == 1
        figure(1);
        plot(abs(autocorr_func));
        xlabel('Starting index');
        ylabel('Autocorrelation amplitude');
        title('Autocorrelation Function');
        grid on;
        set(findall(gcf,'-property','FontSize'),'FontSize',15);
    end
    
    % Fine frequency estimation
    rx_symbol = rx_waveform_delayed(estimated_start + CYCLIC_PREFIX : estimated_start + SUBCARRIERS + CYCLIC_PREFIX - 1);
    fine_freq_offset = angle(sum(rx_symbol(SUBCARRIERS/2+1:end) .* conj(rx_symbol(1:SUBCARRIERS/2)))) / pi;
    
    % Process second synchronization symbol
    rx_symbol_2 = rx_waveform_delayed(estimated_start + CYCLIC_PREFIX*2 + SUBCARRIERS : estimated_start + 2*SUBCARRIERS + 2*CYCLIC_PREFIX - 1);
    rx_symbol_2_corrected = rx_symbol_2 .* exp(-1i * 2 * pi * fine_freq_offset * (0:SUBCARRIERS-1) / SUBCARRIERS);
    rx_symbol_2_fft = fft(rx_symbol_2_corrected);
    
    % Integer frequency estimation
 pn_ratio = tx_symbols(2:2:end, 2) ./ tx_symbols(1:2:end, 2);
q_metric = zeros(1, SUBCARRIERS/2);
for q = 1:SUBCARRIERS/2
    shifted_pn_ratio = circshift(pn_ratio.', q);  % Adjust circshift
    q_metric(q) = abs(sum(rx_symbol_2_fft(2:2:end) .* conj(rx_symbol_2_fft(1:2:end)) .* conj(shifted_pn_ratio))) / sum(abs(rx_symbol_2_fft(2:2:end)).^2);
end
[~, q_est] = max(q_metric);
total_freq_offset = fine_freq_offset + 2 * (q_est - 1);  % Adjust q_est
end

% Plot CDF of timing errors
figure(2);
cdfplot(timing_errors * 1e6);
xlabel('Timing Error (μs)');
ylabel('CDF');
title('CDF of Timing Error');
grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',15);
