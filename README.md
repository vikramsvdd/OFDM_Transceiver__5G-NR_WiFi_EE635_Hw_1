# OFDM_Transceiver__5G-NR_WiFi_EE635_Hw_1

## OFDM Synchronization Simulation
The ofdm_schmidl_cox.m simulates an optimized Schmidl-Cox synchronization algorithm for an OFDM (Orthogonal Frequency-Division Multiplexing) transceiver system. It demonstrates time and frequency synchronization in the presence of noise and carrier frequency offset. The Schiml-Cox model motivated synchronization in 5GNR and WiFi. 
Features

Simulates a 64-subcarrier OFDM system with 1 MHz bandwidth
Implements Schmidl-Cox synchronization algorithm
Handles timing offset and carrier frequency offset
Performs both coarse and fine frequency synchronization
Includes AWGN (Additive White Gaussian Noise) channel simulation
Generates visualizations for autocorrelation function and timing error CDF

## OFDM Transceiver with Walsh Hadamard Transform
This ofdm_hadamard.m simulates a simplified OFDM (Orthogonal Frequency-Division Multiplexing) system using Hadamard transform. It demonstrates the basic principles of OFDM transmission, including modulation, channel effects, and demodulation, while using Hadamard transform instead of just the traditional traditional Fourier transform.
Features

BPSK modulation
Hadamard transform for subcarrier mapping
Cyclic prefix addition and removal
Simulated multipath channel with additive white Gaussian noise (AWGN)
Zero-Forcing equalization
Monte Carlo simulations for Bit Error Rate (BER) calculation
