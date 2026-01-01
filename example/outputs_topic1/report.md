# Topic 1 Report

## Setup
- MATLAB: 25.1.0.2973910 (R2025a) Update 1 on PCWIN64
- Audio: C:\Users\Owner\OneDrive\Homework\DSP\Example\OSR_us_000_0011_8k.wav
- Sample rate: 8000.0 Hz

## EQ Settings
- Crossovers: 300 Hz, 3000 Hz
- Gains: low=+6.0 dB, mid=-3.0 dB, high=+4.0 dB

## Filter Info
- FIR order: 146 (delay: 73 samples = 9.12 ms)
- IIR: Linkwitz-Riley, orders 8/8 at fc1/fc2
- Max pole radius: 0.9087 (stable)

## Reconstruction
- FIR sum error: 1.912e-14
- IIR sum error: 4.247e-15
- Magnitude ripple: FIR=0.000 dB, IIR=0.000 dB

## Echo
- Delay: 250 ms, decay: 0.45, taps: 3

## Output Files
- processed_fir_eq_echo.wav
- processed_iir_eq_echo.wav
- Various PNG plots
