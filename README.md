# DSP Audio Equalizer and Special Effects Processing System

A 3-band crossover equalizer with multi-tap echo effect implemented in MATLAB. Compares FIR (Kaiser window) and IIR (Linkwitz-Riley) filter architectures.

## Features

- 3-band crossover filtering (Low: <300 Hz, Mid: 300-3000 Hz, High: >3000 Hz)
- FIR implementation using Kaiser window method (linear phase, 147 taps)
- IIR implementation using Linkwitz-Riley 4th-order crossovers (flat magnitude reconstruction)
- Multi-tap feedforward echo effect (250ms delay, 3 taps)
- Complete frequency/time-domain analysis and visualization

## Requirements

- MATLAB R2020a or later
- Signal Processing Toolbox

## Usage
```matlab
dsp_final('your_audio_file.wav')
```

Outputs are saved to `outputs_topic1/` directory including:
- Filter magnitude/phase/pole-zero plots
- Spectrograms and waveform comparisons
- Processed audio files (FIR and IIR versions)
- Reconstruction verification plots

## Artifacts

Please avoid committing regenerated files under `outputs_topic1/`; add entries like `outputs_topic1/` (and any similar future output directories) to `.gitignore` when producing new results. The committed examples are provided for demonstration only.

## Technical Details

| Parameter | FIR | IIR |
|-----------|-----|-----|
| Latency | 9.12 ms (73 samples) | Frequency-dependent |
| Phase | Linear | Nonlinear |
| Reconstruction Error | ~10⁻¹⁴ | ~10⁻¹⁵ |

## Example

See the `example/` folder for a complete run using a sample speech file from the Open Speech Repository.

## License

MIT
