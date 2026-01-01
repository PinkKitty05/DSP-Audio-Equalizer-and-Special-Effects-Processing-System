function dsp_final(audioPath)
%% DSP Final Project - 3-Band Equalizer + Echo
% This implements a 3-band crossover EQ using both FIR and IIR approaches
% then adds an echo effect. The IIR version uses Linkwitz-Riley filters.

    %PA RAMETERS 
    % feel free to tweak these
    
    plotDuration = Inf;  % how many seconds to show in plots, this uses the full audio length
    
    % crossover frequencies (Hz)
    p.fc1 = 300;   % splits low from mid
    p.fc2 = 3000;  % splits mid from high
    
    % EQ gains in dB
    p.gain_low_db  = +6;
    p.gain_mid_db  = -3;
    p.gain_high_db = +4;
    
    % FIR settings
    p.fir.Rp_db = 0.5;
    p.fir.Rs_db = 60;         % stopband atten
    p.fir.trans1_hz = 120;    
    p.fir.trans2_hz = 200;    
    
    % IIR Linkwitz-Riley settings
    p.iir.lr_order = 4;       % butterworth order (LR will be 2x this)
    p.iir.Rp_db = 0.5;        
    p.iir.Rs_db = 35;         
    p.iir.guard_hz = 300;     
    p.iir.max_order = 8;      
    p.iir.min_Rs_db = 30;     
    p.iir.max_guard_hz = 600; 
    
    % echo settings
    p.echo.delay_ms = 250;
    p.echo.alpha = 0.45;   % decay factor
    p.echo.numTaps = 3;    % number of echoes
    
    % spectrogram params
    p.spec.win_ms = 30;
    p.spec.overlap = 0.75;
    p.spec.nfft = 4096;
    
    p.playback = false;  % set true to hear output
    
    % SETUP
    
    % make output folder
    outDir = fullfile(pwd, "outputs_topic1");
    if ~exist(outDir, "dir"), mkdir(outDir); end
    clearOldOutputs(outDir);
    
    % get audio file
    if nargin < 1 || strlength(string(audioPath)) == 0
        [fName, pName] = uigetfile({'*.wav;*.flac','Audio'; '*.*','All'}, 'Pick audio file');
        if isequal(fName, 0), error("no file selected"); end
        audioPath = fullfile(pName, fName);
    else
        audioPath = char(audioPath);
    end
    
    if exist(audioPath, "file") ~= 2
        error("cant find file: %s", audioPath);
    end
    
    % load audio
    [x, fs] = audioread(audioPath);
    x = makeMono(x);
    x = x / max(abs(x));  % normalize
    
    % for plotting just use first few seconds
    Nfull = numel(x);
    Nplot = min(Nfull, round(plotDuration * fs));
    xPlot = x(1:Nplot);
    xPlot = xPlot - mean(xPlot);  % remove DC
    tPlot = (0:Nplot-1) / fs;
    
    % make sure crossovers make sense
    if ~(p.fc1 > 0 && p.fc1 < p.fc2 && p.fc2 < fs/2)
        error("crossover freqs dont work with this sample rate");
    end
    
    % store some info for the summary
    info.audioPath = string(audioPath);
    info.fs = fs;
    info.nyquist = fs/2;
    info.matlabVersion = version;
    info.platform = computer;
    
    % PLOT OG
    
    fig = figure('Visible', 'off');
    subplot(2,1,1)
    plot(tPlot, xPlot); grid on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Original waveform');
    ylim([-1 1]);
    
    subplot(2,1,2)
    % zoom in on a small chunk
    t0 = 2.0; dur = 0.05; 
    i0 = max(1, round(t0*fs));
    i1 = min(numel(xPlot), i0 + round(dur*fs));
    plot(tPlot(i0:i1), xPlot(i0:i1)); grid on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Zoomed (50 ms around t=2s)');
    ylim([-1 1]);
    savePlot(fig, outDir, "01_original_waveform.png");
    
    fig = figure('Visible', 'off');
    doSpectrum(xPlot, fs);
    title('Original spectrum');
    savePlot(fig, outDir, "02_original_spectrum.png");
    
    fig = figure('Visible', 'off');
    doSpectrogram(xPlot, fs, p.spec);
    title('Original spectrogram');
    savePlot(fig, outDir, "03_original_spectrogram.png");
    
    % DESIGN FIR FILTERS
    
    [bF_L, bF_M, bF_H, firOrder, firDelay] = makeFIR3Band(fs, p.fc1, p.fc2, p.fir);
    
    % convert dB gains to linear
    gL = 10^(p.gain_low_db/20);
    gM = 10^(p.gain_mid_db/20);
    gH = 10^(p.gain_high_db/20);
    
    % filter and combine
    yF_L = filter(bF_L, 1, x);
    yF_M = filter(bF_M, 1, x);
    yF_H = filter(bF_H, 1, x);
    yEQ_FIR = gL*yF_L + gM*yF_M + gH*yF_H;
    
    % check reconstruction (should sum to delayed impulse)
    recon.nfft = max(2048, p.spec.nfft);
    [recon.fir_f_hz, recon.fir_Hsum, recon.fir_Htarget] = checkFIRReconstruction(bF_L, bF_M, bF_H, firDelay, fs, recon.nfft);
    
    % run some sanity checks
    firChecks = struct();
    [firChecks.low_dc, firChecks.low_nyq] = getEdgeGains(bF_L);
    [firChecks.mid_dc, firChecks.mid_nyq] = getEdgeGains(bF_M);
    [firChecks.high_dc, firChecks.high_nyq] = getEdgeGains(bF_H);
    [firChecks.low_sym_abs, firChecks.low_sym_rel] = checkSymmetry(bF_L);
    [firChecks.mid_sym_abs, firChecks.mid_sym_rel] = checkSymmetry(bF_M);
    [firChecks.high_sym_abs, firChecks.high_sym_rel] = checkSymmetry(bF_H);
    
    % group delay stuff
    wF = 2*pi*(recon.fir_f_hz(:)/fs);
    phiF = unwrap(angle(recon.fir_Hsum(:)));
    gdF = -diff(phiF) ./ (diff(wF) + 1e-12);
    firChecks.hsum_gd_median = median(gdF(~isnan(gdF) & ~isinf(gdF)));
    firChecks.hsum_gd_max_abs_dev = max(abs(gdF - firChecks.hsum_gd_median));
    firChecks.hsum_mag_max_abs_dev_db = max(abs(20*log10(abs(recon.fir_Hsum(:)) + 1e-12)));
    
    [firChecks.hsum_best_delay, firChecks.hsum_best_sign, firChecks.hsum_best_max_err, firChecks.hsum_best_rms_err] = ...
        findBestDelay(recon.fir_Hsum, wF, max(50, round(2*firDelay + 10)));
    
    info.fir_checks = firChecks;
    info.fir_order = firOrder;
    info.fir_group_delay_samples = firDelay;
    info.fir_latency_ms = 1000*firDelay/fs;
    
    % FIR PLOTS
    
    fig = figure('Visible', 'off');
    plotMagFIR(bF_L, 1, fs); title("FIR Low: magnitude");
    savePlot(fig, outDir, "04_fir_01_low_magnitude.png");
    
    fig = figure('Visible', 'off');
    plotPhaseFIR(bF_L, 1, fs); title("FIR Low: phase");
    savePlot(fig, outDir, "04_fir_02_low_phase.png");
    
    fig = figure('Visible', 'off');
    plotZplane(gca, bF_L, 1); title("FIR Low: z plane");
    savePlot(fig, outDir, "04_fir_03_low_zplane.png");
    
    fig = figure('Visible', 'off');
    plotMagFIR(bF_M, 1, fs); title("FIR Mid: magnitude");
    savePlot(fig, outDir, "04_fir_04_mid_magnitude.png");
    
    fig = figure('Visible', 'off');
    plotPhaseFIR(bF_M, 1, fs); title("FIR Mid: phase");
    savePlot(fig, outDir, "04_fir_05_mid_phase.png");
    
    fig = figure('Visible', 'off');
    plotZplane(gca, bF_M, 1); title("FIR Mid: z plane");
    savePlot(fig, outDir, "04_fir_06_mid_zplane.png");
    
    fig = figure('Visible', 'off');
    plotMagFIR(bF_H, 1, fs); title("FIR High: magnitude");
    savePlot(fig, outDir, "04_fir_07_high_magnitude.png");
    
    fig = figure('Visible', 'off');
    plotPhaseFIR(bF_H, 1, fs); title("FIR High: phase");
    savePlot(fig, outDir, "04_fir_08_high_phase.png");
    
    fig = figure('Visible', 'off');
    plotZplane(gca, bF_H, 1); title("FIR High: z plane");
    savePlot(fig, outDir, "04_fir_09_high_zplane.png");
    
    % DESIGN IIR (LINKWITZ-RILEY) 
    
    [sL1, gL1, sH1, gH1, sL2, gL2, sH2, gH2, iirOrders, iirUsed] = makeIIR3BandLR(fs, p.fc1, p.fc2, p.iir);
    
    % The 3-band LR crossover works like this:
    % - first split at fc1 into low1 and high1
    % - then split high1 at fc2 into mid and high
    % - for phase alignment, low band needs to go through the second crossover too
    
    y_L1 = applySOS(sL1, gL1, x);
    y_H1 = applySOS(sH1, gH1, x);
    
    % low band gets phase compensation by going thru crossover 2
    y_L1_L2 = applySOS(sL2, gL2, y_L1);
    y_L1_H2 = applySOS(sH2, gH2, y_L1);
    yI_L = y_L1_L2 + y_L1_H2;  % this is the phase-aligned low band
    
    % mid = HP1 * LP2, high = HP1 * HP2
    yI_M = applySOS(sL2, gL2, y_H1);
    yI_H = applySOS(sH2, gH2, y_H1);
    
    % combine with gains
    yEQ_IIR = gL*yI_L + gM*yI_M + gH*yI_H;
    
    % for plotting purposes
    sosI_L = sL1;           gI_L = gL1;
    sosI_M = [sH1; sL2];    gI_M = gH1 * gL2;
    sosI_H = [sH1; sH2];    gI_H = gH1 * gH2;
    
    info.iir_orders = [iirOrders(1) iirOrders(3)];
    info.iir_Rp_db_used = iirUsed.Rp_db;
    info.iir_Rs_db_used = iirUsed.Rs_db;
    info.iir_guard_hz_used = iirUsed.guard_hz;
    
    % stability check
    [~, pL] = getZP(sosI_L, gI_L);
    [~, pM] = getZP(sosI_M, gI_M);
    [~, pH] = getZP(sosI_H, gI_H);
    allPoles = [pL(:); pM(:); pH(:)];
    info.iir_max_pole_radius = max(abs(allPoles));
    
    % quick sanity check - LR crossovers should sum to unity
    [HL1_test, ~] = getFreqResp(sL1, gL1, 2048);
    HH1_test = getFreqResp(sH1, gH1, 2048);
    sum1 = abs(HL1_test + HH1_test);
    fprintf('=== LR Crossover Checks ===\n');
    fprintf('fc1=%d Hz: sum ranges from %.6f to %.6f (should be ~1)\n', p.fc1, min(sum1), max(sum1));
    
    [HL2_test, ~] = getFreqResp(sL2, gL2, 2048);
    HH2_test = getFreqResp(sH2, gH2, 2048);
    sum2 = abs(HL2_test + HH2_test);
    fprintf('fc2=%d Hz: sum ranges from %.6f to %.6f (should be ~1)\n', p.fc2, min(sum2), max(sum2));
    
    % IIR RECONSTRUCTION CHECK
    
    [HL1, wI] = getFreqResp(sL1, gL1, recon.nfft);
    HH1 = getFreqResp(sH1, gH1, recon.nfft);
    HL2 = getFreqResp(sL2, gL2, recon.nfft);
    HH2 = getFreqResp(sH2, gH2, recon.nfft);
    
    % the actual band transfer functions
    H_L_comp = HL1 .* (HL2 + HH2);
    H_M_comp = HH1 .* HL2;
    H_H_comp = HH1 .* HH2;
    
    recon.iir_Hsum = H_L_comp + H_M_comp + H_H_comp;
    recon.iir_f_hz = (wI/(2*pi))*fs;
    
    % calculate errors
    recon.fir_max_mag_err_db = max(20*log10(abs(recon.fir_Hsum) + 1e-12));
    recon.fir_min_mag_err_db = min(20*log10(abs(recon.fir_Hsum) + 1e-12));
    recon.fir_max_abs_dev = max(abs(recon.fir_Hsum - recon.fir_Htarget));
    
    recon.iir_max_mag_err_db = max(20*log10(abs(recon.iir_Hsum) + 1e-12));
    recon.iir_min_mag_err_db = min(20*log10(abs(recon.iir_Hsum) + 1e-12));
    recon.iir_max_abs_dev = max(abs(recon.iir_Hsum - 1));
    
    magF = abs(recon.fir_Hsum(:));
    recon.fir_mag_min = min(magF);
    recon.fir_mag_max = max(magF);
    recon.fir_mag_ripple_db = 20*log10(recon.fir_mag_max + 1e-12) - 20*log10(recon.fir_mag_min + 1e-12);
    
    magI = abs(recon.iir_Hsum(:));
    recon.iir_mag_min = min(magI);
    recon.iir_mag_max = max(magI);
    recon.iir_mag_ripple_db = 20*log10(recon.iir_mag_max + 1e-12) - 20*log10(recon.iir_mag_min + 1e-12);
    
    % time domain check
    y_ap2 = applySOS(sL2, gL2, x) + applySOS(sH2, gH2, x);
    y_target = applySOS(sL1, gL1, y_ap2) + applySOS(sH1, gH1, y_ap2);
    recon.iir_time_recon_max_err = max(abs((yI_L + yI_M + yI_H) - y_target));
    
    % alternative check using residual
    yI_M_resid = x - yI_L - yI_H;
    recon.iir_time_recon_max_err_resid = max(abs((yI_L + yI_M_resid + yI_H) - x));
    
    info.recon = recon;
    
    % IIR PLOTS 
    % note: the phase responses end up identical across bands
    % this is actually correct for LR - its how they sum flat
    
    H_Low  = HL1 .* (HL2 + HH2);
    H_Mid  = HH1 .* HL2;
    H_High = HH1 .* HH2;
    
    % low band
    fig = figure('Visible', 'off');
    plot(recon.iir_f_hz, 20*log10(abs(H_Low) + 1e-12), 'LineWidth', 1.2); grid on;
    xlim([0 fs/2]); ylim([-80 10]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title("IIR Low Band: magnitude");
    savePlot(fig, outDir, "05_iir_01_low_magnitude.png");
    
    fig = figure('Visible', 'off');
    plot(recon.iir_f_hz, unwrap(angle(H_Low))*180/pi, 'LineWidth', 1.2); grid on;
    xlim([0 fs/2]); xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
    title("IIR Low Band: phase");
    savePlot(fig, outDir, "05_iir_02_low_phase.png");
    
    fig = figure('Visible', 'off');
    plotZplaneSOS(gca, sL1, gL1);
    title(sprintf("IIR Low Band (LP @ %d Hz): z-plane", p.fc1));
    savePlot(fig, outDir, "05_iir_03_low_zplane.png");
    
    % mid band
    fig = figure('Visible', 'off');
    plot(recon.iir_f_hz, 20*log10(abs(H_Mid) + 1e-12), 'LineWidth', 1.2); grid on;
    xlim([0 fs/2]); ylim([-80 10]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title("IIR Mid Band: magnitude");
    savePlot(fig, outDir, "05_iir_04_mid_magnitude.png");
    
    fig = figure('Visible', 'off');
    plot(recon.iir_f_hz, unwrap(angle(H_Mid))*180/pi, 'LineWidth', 1.2); grid on;
    xlim([0 fs/2]); xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
    title("IIR Mid Band: phase");
    savePlot(fig, outDir, "05_iir_05_mid_phase.png");
    
    fig = figure('Visible', 'off');
    plotZplaneSOS(gca, [sH1; sL2], gH1*gL2);
    title("IIR Mid Band (HP1*LP2): z-plane");
    savePlot(fig, outDir, "05_iir_06_mid_zplane.png");
    
    % high band
    fig = figure('Visible', 'off');
    plot(recon.iir_f_hz, 20*log10(abs(H_High) + 1e-12), 'LineWidth', 1.2); grid on;
    xlim([0 fs/2]); ylim([-80 10]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title("IIR High Band: magnitude");
    savePlot(fig, outDir, "05_iir_07_high_magnitude.png");
    
    fig = figure('Visible', 'off');
    plot(recon.iir_f_hz, unwrap(angle(H_High))*180/pi, 'LineWidth', 1.2); grid on;
    xlim([0 fs/2]); xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
    title("IIR High Band: phase");
    savePlot(fig, outDir, "05_iir_08_high_phase.png");
    
    fig = figure('Visible', 'off');
    plotZplaneSOS(gca, [sH1; sH2], gH1*gH2);
    title("IIR High Band (HP1*HP2): z-plane");
    savePlot(fig, outDir, "05_iir_09_high_zplane.png");
    
    % reconstruction plots
    fig = figure('Visible', 'off');
    plot(recon.fir_f_hz, 20*log10(abs(recon.fir_Hsum) + 1e-12)); grid on;
    xlabel('Frequency (Hz)'); ylabel('Magnitude |H_L+H_M+H_H| (dB)');
    title('FIR reconstruction (should be flat at 0 dB)');
    xlim([0 fs/2]);
    savePlot(fig, outDir, "10_reconstruction_01_fir_hsum.png");
    
    fig = figure('Visible', 'off');
    plot(recon.iir_f_hz, 20*log10(abs(recon.iir_Hsum) + 1e-12)); grid on;
    xlabel('Frequency (Hz)'); ylabel('Magnitude |H_L+H_M+H_H| (dB)');
    title('IIR LR reconstruction (should be flat at 0 dB)');
    xlim([0 fs/2]);
    savePlot(fig, outDir, "10_reconstruction_02_iir_hsum.png");
    
    % ECHO + HEADROOM
    
    yEcho_FIR = addEcho(yEQ_FIR, fs, p.echo.delay_ms, p.echo.alpha, p.echo.numTaps);
    yEcho_IIR = addEcho(yEQ_IIR, fs, p.echo.delay_ms, p.echo.alpha, p.echo.numTaps);
    
    % prevent clipping
    [yEcho_FIR, info.fir_peak_before, info.fir_peak_after] = fixHeadroom(yEcho_FIR, -1);
    [yEcho_IIR, info.iir_peak_before, info.iir_peak_after] = fixHeadroom(yEcho_IIR, -1);
    
    % match loudness to original
    yEcho_FIR = matchLoudness(yEcho_FIR, x);
    yEcho_IIR = matchLoudness(yEcho_IIR, x);
    
    info.fir_peak_final = max(abs(yEcho_FIR));
    info.iir_peak_final = max(abs(yEcho_IIR));
    
    % stats
    stats.rms_x = rms(x);
    stats.rms_fir = rms(yEcho_FIR);
    stats.rms_iir = rms(yEcho_IIR);
    stats.crest_x = max(abs(x)) / (rms(x) + 1e-12);
    stats.crest_fir = max(abs(yEcho_FIR)) / (rms(yEcho_FIR) + 1e-12);
    stats.crest_iir = max(abs(yEcho_IIR)) / (rms(yEcho_IIR) + 1e-12);
    
    % COMPARISON PLOT
    % need to shift FIR output to compensate for group delay
    yPlot_FIR = shiftForPlot(yEcho_FIR(1:Nplot), firDelay);
    yPlot_IIR = yEcho_IIR(1:Nplot);
    
    % time domain
    fig = figure('Visible', 'off');
    plot(tPlot, yPlot_FIR); grid on; title('FIR EQ + Echo');
    xlabel('Time (s)'); ylabel('Amplitude');
    savePlot(fig, outDir, "06_waveform_02_fir_eq_echo.png");
    
    fig = figure('Visible', 'off');
    plot(tPlot, yPlot_IIR); grid on; title('IIR EQ + Echo');
    xlabel('Time (s)'); ylabel('Amplitude');
    savePlot(fig, outDir, "06_waveform_03_iir_eq_echo.png");
    
    % spectra
    specCfg = struct('nfft', p.spec.nfft);
    refSpec = getRefSpectrum(xPlot, fs, specCfg);
    
    fig = figure('Visible', 'off');
    doSpectrum(yPlot_FIR, fs, specCfg, refSpec); title('FIR EQ + Echo spectrum');
    savePlot(fig, outDir, "07_spectrum_01_fir_eq_echo.png");
    
    fig = figure('Visible', 'off');
    doSpectrum(yPlot_IIR, fs, specCfg, refSpec); title('IIR EQ + Echo spectrum');
    savePlot(fig, outDir, "07_spectrum_02_iir_eq_echo.png");
    
    % spectrograms
    fig = figure('Visible', 'off');
    doSpectrogram(yPlot_FIR, fs, p.spec); title('FIR EQ + Echo spectrogram');
    savePlot(fig, outDir, "08_spectrogram_01_fir_eq_echo.png");
    
    fig = figure('Visible', 'off');
    doSpectrogram(yPlot_IIR, fs, p.spec); title('IIR EQ + Echo spectrogram');
    savePlot(fig, outDir, "08_spectrogram_02_iir_eq_echo.png");
    
    % echo filter response
    fig = figure('Visible', 'off');
    plotEchoMag(fs, p.echo.delay_ms, p.echo.alpha, p.echo.numTaps);
    savePlot(fig, outDir, "09_echo_01_magnitude.png");
    
    fig = figure('Visible', 'off');
    plotEchoPhase(fs, p.echo.delay_ms, p.echo.alpha, p.echo.numTaps);
    savePlot(fig, outDir, "09_echo_02_phase.png");
    
    % SAVE OUTPUTS
    
    audiowrite(fullfile(outDir, "processed_fir_eq_echo.wav"), yEcho_FIR, fs);
    audiowrite(fullfile(outDir, "processed_iir_eq_echo.wav"), yEcho_IIR, fs);
    
    % optional playback
    if p.playback && usejava('desktop')
        try
            disp('Playing original...');
            soundsc(xPlot, fs); pause(numel(xPlot)/fs + 0.25);
            disp('Playing FIR...');
            soundsc(yPlot_FIR, fs); pause(numel(yPlot_FIR)/fs + 0.25);
            disp('Playing IIR...');
            soundsc(yPlot_IIR, fs); pause(numel(yPlot_IIR)/fs + 0.25);
        catch
        end
    end
    
    % write summary
    writeSummary(fullfile(outDir, "run_summary.txt"), info, stats, p);
    
    disp("Done! Check " + string(outDir));
end

function x = makeMono(x)
    if size(x,2) > 1
        x = mean(x, 2);
    end
    x = x(:);
end

function [bL, bM, bH, N, gd] = makeFIR3Band(fs, fc1, fc2, spec)
    % design FIR crossover using kaiser window method
    Rs = spec.Rs_db;
    tw = max([spec.trans1_hz spec.trans2_hz]);
    
    A = Rs;
    if A > 50
        beta = 0.1102*(A - 8.7);
    elseif A >= 21
        beta = 0.5842*(A - 21)^0.4 + 0.07886*(A - 21);
    else
        beta = 0;
    end
    
    dw = 2*pi*(tw/fs);
    N = ceil((A - 8) / (2.285*dw));
    if mod(N,2)==1, N = N + 1; end
    N = max(N, 64);
    
    win = kaiser(N+1, beta);
    
    bL = fir1(N, fc1/(fs/2), 'low', win, 'scale');
    bH = fir1(N, fc2/(fs/2), 'high', win, 'scale');
    gd = N/2;
    
    % mid is the complement so they sum to delayed impulse
    delta = zeros(1, N+1);
    delta(gd+1) = 1;
    bM = delta - bL - bH;
end

function [sLP1, gLP1, sHP1, gHP1, sLP2, gLP2, sHP2, gHP2, orders, used] = makeIIR3BandLR(fs, fc1, fc2, spec)
    % design linkwitz-riley 3-band crossover
    Rp = spec.Rp_db;
    Rs = spec.Rs_db;
    g  = spec.guard_hz;
    maxN = spec.max_order;
    minRs = spec.min_Rs_db;
    maxG = spec.max_guard_hz;
    
    wp1 = fc1 / (fs/2);
    wp2 = fc2 / (fs/2);
    
    % if lr_order specified, use that directly
    if isfield(spec, 'lr_order') && ~isempty(spec.lr_order)
        lrN = max(1, round(spec.lr_order));
        
        % tune cutoff for best flatness
        grid = linspace(0.92, 1.08, 61);
        nTune = 4096;
        
        [sLP1, gLP1, sHP1, gHP1] = tuneLRPair(lrN, wp1, grid, nTune);
        [sLP2, gLP2, sHP2, gHP2] = tuneLRPair(lrN, wp2, grid, nTune);
        
        orders = [2*lrN 2*lrN 2*lrN 2*lrN];
        used.Rp_db = Rp; used.Rs_db = Rs; used.guard_hz = g; used.lr_order = lrN;
        return;
    end
    
    % otherwise figure out orders automatically
    % first crossover
    n1 = NaN; Wn1 = wp1;
    Rs1 = Rs; g1 = g;
    for k = 1:12
        wsLP = min((fc1 + g1)/(fs/2), 0.99);
        wsHP = max((fc1 - g1)/(fs/2), 0.01);
        [nLP, ~] = buttord(wp1, wsLP, Rp, Rs1);
        [nHP, ~] = buttord(wp1, wsHP, Rp, Rs1);
        n1 = max(nLP, nHP);
        if mod(n1, 2) == 1, n1 = n1 + 1; end
        if n1 <= maxN, break; end
        if g1 < maxG
            g1 = min(maxG, round(g1 * 1.4));
        elseif Rs1 > minRs
            Rs1 = max(minRs, Rs1 - 5);
        else
            n1 = maxN; if mod(n1,2)==1, n1=n1-1; end
            break;
        end
    end
    n1 = max(2, min(n1, maxN));
    if mod(n1,2)==1, n1 = n1 + (n1<maxN) - (n1>=maxN); end
    Wn1 = min(max(wp1, 1e-4), 0.999);
    
    % second crossover (same process)
    n2 = NaN; Wn2 = wp2;
    Rs2 = Rs; g2 = g;
    for k = 1:12
        wsLP = min((fc2 + g2)/(fs/2), 0.99);
        wsHP = max((fc2 - g2)/(fs/2), 0.01);
        [nLP, ~] = buttord(wp2, wsLP, Rp, Rs2);
        [nHP, ~] = buttord(wp2, wsHP, Rp, Rs2);
        n2 = max(nLP, nHP);
        if mod(n2, 2) == 1, n2 = n2 + 1; end
        if n2 <= maxN, break; end
        if g2 < maxG
            g2 = min(maxG, round(g2 * 1.4));
        elseif Rs2 > minRs
            Rs2 = max(minRs, Rs2 - 5);
        else
            n2 = maxN; if mod(n2,2)==1, n2=n2-1; end
            break;
        end
    end
    n2 = max(2, min(n2, maxN));
    if mod(n2,2)==1, n2 = n2 + (n2<maxN) - (n2>=maxN); end
    Wn2 = min(max(wp2, 1e-4), 0.999);
    
    % design butterworth then cascade for LR
    [sLP1_b, gLP1_b] = getButter(n1, Wn1, 'low');
    [sHP1_b, gHP1_b] = getButter(n1, Wn1, 'high');
    [sLP2_b, gLP2_b] = getButter(n2, Wn2, 'low');
    [sHP2_b, gHP2_b] = getButter(n2, Wn2, 'high');
    
    sLP1 = [sLP1_b; sLP1_b]; gLP1 = gLP1_b^2;
    sHP1 = [sHP1_b; sHP1_b]; gHP1 = gHP1_b^2;
    sLP2 = [sLP2_b; sLP2_b]; gLP2 = gLP2_b^2;
    sHP2 = [sHP2_b; sHP2_b]; gHP2 = gHP2_b^2;
    
    orders = [n1 n1 n2 n2];
    used.Rp_db = Rp;
    used.Rs_db = max(Rs1, Rs2);
    used.guard_hz = max(g1, g2);
end

function [sLP, gLP, sHP, gHP] = tuneLRPair(n, wp, grid, nfft)
    % find best cutoff frequency for flat LR response
    wp = min(max(wp, 1e-6), 0.999);
    bestErr = Inf;
    best = struct('sLP',[],'gLP',[],'sHP',[],'gHP',[]);
    
    for k = 1:numel(grid)
        Wn = min(max(wp * grid(k), 1e-6), 0.999);
        [sLP_b, gLP_b] = getButter(n, Wn, 'low');
        [sHP_b, gHP_b] = getButter(n, Wn, 'high');
        
        sLP_k = [sLP_b; sLP_b]; gLP_k = gLP_b^2;
        sHP_k = [sHP_b; sHP_b]; gHP_k = gHP_b^2;
        
        HL = getFreqResp(sLP_k, gLP_k, nfft);
        HH = getFreqResp(sHP_k, gHP_k, nfft);
        err = max(abs(abs(HL + HH) - 1));
        
        if err < bestErr
            bestErr = err;
            best.sLP = sLP_k; best.gLP = gLP_k;
            best.sHP = sHP_k; best.gHP = gHP_k;
        end
    end
    sLP = best.sLP; gLP = best.gLP;
    sHP = best.sHP; gHP = best.gHP;
end

function y = applySOS(sos, g, x)
    x = x(:);
    if exist('sosfilt', 'file') == 2
        y = g * sosfilt(sos, x);
    else
        y = x;
        for i = 1:size(sos,1)
            y = filter(sos(i,1:3), sos(i,4:6), y);
        end
        y = g * y;
    end
end

function [sos, g] = getButter(n, Wn, ftype)
    try
        [sos, g] = butter(n, Wn, ftype, 'sos');
    catch
        [z, p, k] = butter(n, Wn, ftype);
        [sos, g] = zp2sos(z, p, k);
    end
end

function [z, p] = getZP(sos, g)
    if exist('sos2zp', 'file') == 2
        [z, p] = sos2zp(sos, g);
    else
        z = []; p = [];
        for i = 1:size(sos,1)
            z = [z; roots(sos(i,1:3))];
            p = [p; roots(sos(i,4:6))];
        end
    end
end

function [H, w] = getFreqResp(sos, g, nfft)
    if nargin < 3, nfft = 2048; end
    
    if exist('sosfreqz', 'file') == 2
        [H, w] = sosfreqz(sos, nfft);
        H = g * H;
    else
        w = linspace(0, pi, nfft);
        H = ones(size(w));
        for i = 1:size(sos,1)
            H = H .* freqz(sos(i,1:3), sos(i,4:6), w);
        end
        H = g * H;
    end
end

function [dcG, nyqG] = getEdgeGains(b)
    b = b(:);
    n = (0:numel(b)-1).';
    dcG = abs(sum(b));
    nyqG = abs(sum(b .* ((-1).^n)));
end

function [symAbs, symRel] = checkSymmetry(b)
    b = b(:);
    symAbs = max(abs(b - flipud(b)));
    symRel = symAbs / (max(abs(b)) + 1e-12);
end

function [f, Hsum, Htarget] = checkFIRReconstruction(bL, bM, bH, gd, fs, nfft)
    [HL, w] = freqz(bL, 1, nfft);
    HM = freqz(bM, 1, w);
    HH = freqz(bH, 1, w);
    Hsum = HL + HM + HH;
    Htarget = exp(-1j*w*gd);
    f = (w/(2*pi))*fs;
end

function [bestD, bestSign, maxErr, rmsErr] = findBestDelay(H, w, dMax)
    H = H(:); w = w(:);
    if isempty(H), bestD=0; bestSign=1; maxErr=NaN; rmsErr=NaN; return; end
    
    bestD = 0; bestSign = 1; bestRms = inf; maxErr = inf;
    for s = [1 -1]
        for d = 0:dMax
            T = s * exp(-1j*w*d);
            e = H - T;
            r = sqrt(mean(abs(e).^2));
            if r < bestRms
                bestRms = r; bestD = d; bestSign = s;
                maxErr = max(abs(e));
            end
        end
    end
    rmsErr = bestRms;
end

function y = addEcho(x, fs, delay_ms, alpha, numTaps)
    x = x(:);
    delaySamp = round(delay_ms/1000 * fs);
    
    b = zeros(1, numTaps*delaySamp + 1);
    for k = 0:numTaps-1
        b(k*delaySamp + 1) = alpha^k;
    end
    b = b / sum(abs(b));
    
    y = filter(b, 1, x);
end

function [y, peakBefore, peakAfter] = fixHeadroom(x, targetDb)
    x = x(:);
    peakBefore = max(abs(x));
    if targetDb < 0
        targetLin = 10^(targetDb/20);
    else
        targetLin = targetDb;
    end
    if peakBefore > targetLin
        y = x * (targetLin / peakBefore);
    else
        y = x;
    end
    peakAfter = max(abs(y));
end

function y = matchLoudness(x, ref)
    rmsRef = rms(ref);
    rmsX = rms(x);
    if rmsX > 1e-12
        y = x * (rmsRef / rmsX);
    else
        y = x;
    end
end

function y = shiftForPlot(x, d)
    x = x(:);
    N = numel(x);
    if d > 0 && d < N
        y = [x(d+1:end); zeros(d, 1)];
    else
        y = x;
    end
end

function plotMagFIR(b, a, fs)
    [H, w] = freqz(b, a, 2048);
    f = (w/(2*pi))*fs;
    plot(f, 20*log10(abs(H)+1e-12)); grid on;
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); xlim([0 fs/2]);
end

function plotPhaseFIR(b, a, fs)
    [H, w] = freqz(b, a, 2048);
    f = (w/(2*pi))*fs;
    plot(f, unwrap(angle(H))*180/pi); grid on;
    xlabel('Frequency (Hz)'); ylabel('Phase (degrees)'); xlim([0 fs/2]);
end

function plotZplane(ax, b, a)
    z = roots(b(:)); p = roots(a(:));
    cla(ax); hold(ax, 'on');
    th = linspace(0, 2*pi, 400);
    plot(ax, cos(th), sin(th), '--');
    if ~isempty(z), plot(ax, real(z), imag(z), 'o', 'MarkerSize', 6); end
    if ~isempty(p), plot(ax, real(p), imag(p), 'x', 'MarkerSize', 6); end
    axis(ax, 'equal'); grid(ax, 'on');
    xlabel(ax, 'Real'); ylabel(ax, 'Imaginary');
    hold(ax, 'off');
end

function plotZplaneSOS(ax, sos, g)
    [z, p] = getZP(sos, g);
    cla(ax); hold(ax, 'on');
    th = linspace(0, 2*pi, 400);
    plot(ax, cos(th), sin(th), '--');
    if ~isempty(z), plot(ax, real(z), imag(z), 'o', 'MarkerSize', 6); end
    if ~isempty(p), plot(ax, real(p), imag(p), 'x', 'MarkerSize', 6); end
    axis(ax, 'equal'); grid(ax, 'on');
    xlabel(ax, 'Real'); ylabel(ax, 'Imaginary');
    hold(ax, 'off');
end

function doSpectrum(x, fs, cfg, ref)
    if nargin < 3 || isempty(cfg), cfg = struct('nfft', 4096); end
    nfft = cfg.nfft;
    x = x(:);
    X = fft(x .* hann(numel(x)), nfft);
    f = (0:nfft/2-1) * fs / nfft;
    magDb = 20*log10(abs(X(1:nfft/2)) + 1e-12);
    if nargin >= 4 && ~isempty(ref)
        magDb = magDb - ref;
    end
    plot(f, magDb); grid on;
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); xlim([0 fs/2]);
end

function ref = getRefSpectrum(x, fs, cfg)
    nfft = cfg.nfft;
    x = x(:);
    X = fft(x .* hann(numel(x)), nfft);
    ref = 20*log10(abs(X(1:nfft/2)) + 1e-12);
end

function doSpectrogram(x, fs, cfg)
    winLen = round(cfg.win_ms/1000 * fs);
    noverlap = round(cfg.overlap * winLen);
    nfft = cfg.nfft;
    spectrogram(x, hann(winLen), noverlap, nfft, fs, 'yaxis');
    colorbar;
end

function plotEchoMag(fs, delay_ms, alpha, numTaps)
    delaySamp = round(delay_ms/1000 * fs);
    b = zeros(1, numTaps*delaySamp + 1);
    for k = 0:numTaps-1
        b(k*delaySamp + 1) = alpha^k;
    end
    b = b / sum(abs(b));
    [H, w] = freqz(b, 1, 2048);
    f = (w/(2*pi))*fs;
    plot(f, 20*log10(abs(H)+1e-12)); grid on;
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title('Echo filter magnitude');
    xlim([0 fs/2]);
end

function plotEchoPhase(fs, delay_ms, alpha, numTaps)
    delaySamp = round(delay_ms/1000 * fs);
    b = zeros(1, numTaps*delaySamp + 1);
    for k = 0:numTaps-1
        b(k*delaySamp + 1) = alpha^k;
    end
    b = b / sum(abs(b));
    [H, w] = freqz(b, 1, 2048);
    f = (w/(2*pi))*fs;
    plot(f, unwrap(angle(H))*180/pi); grid on;
    xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
    title('Echo filter phase');
    xlim([0 fs/2]);
end

function savePlot(fig, folder, name)
    fpath = fullfile(folder, name);
    try
        exportgraphics(fig, fpath, 'Resolution', 200);
    catch
        print(fig, fpath, '-dpng', '-r200');
    end
    close(fig);
end

function clearOldOutputs(folder)
    patterns = {'*.png', '*.wav', '*.txt', 'report.md'};
    for i = 1:numel(patterns)
        files = dir(fullfile(folder, patterns{i}));
        for k = 1:numel(files)
            try delete(fullfile(folder, files(k).name)); catch, end
        end
    end
end

function writeSummary(path, info, stats, p)
    fid = fopen(path, 'w');
    if fid < 0, return; end
    
    fprintf(fid, 'Topic 1 summary\n\n');
    fprintf(fid, 'Audio: %s\n', info.audioPath);
    fprintf(fid, 'fs: %.1f Hz (Nyquist %.1f Hz)\n', info.fs, info.nyquist);
    fprintf(fid, 'MATLAB: %s on %s\n\n', info.matlabVersion, info.platform);
    
    fprintf(fid, 'Crossovers: fc1=%d Hz, fc2=%d Hz\n', p.fc1, p.fc2);
    fprintf(fid, 'Gains (dB): low=%+.1f, mid=%+.1f, high=%+.1f\n\n', p.gain_low_db, p.gain_mid_db, p.gain_high_db);
    fprintf(fid, 'Preprocessing: mono (channel average when needed), peak normalized\n\n');
    
    fprintf(fid, 'FIR: order=%d, group delay=%d samples (%.2f ms)\n', info.fir_order, info.fir_group_delay_samples, info.fir_latency_ms);
    fprintf(fid, 'IIR (3-band Linkwitz-Riley): base orders @fc1/@fc2=%d/%d\n', info.iir_orders(1), info.iir_orders(2));
    
    if isfield(p.iir, 'lr_order') && ~isempty(p.iir.lr_order)
        fprintf(fid, 'IIR mode: explicit LR (lr_order=%d); Rp/Rs/guard fields are not used\n', round(p.iir.lr_order));
    end
    fprintf(fid, 'IIR max pole radius=%.4f\n\n', info.iir_max_pole_radius);
    
    fprintf(fid, 'Peaks (abs): FIR before headroom=%.4f, after headroom=%.4f, final=%.4f\n', info.fir_peak_before, info.fir_peak_after, info.fir_peak_final);
    fprintf(fid, 'Peaks (abs): IIR before headroom=%.4f, after headroom=%.4f, final=%.4f\n\n', info.iir_peak_before, info.iir_peak_after, info.iir_peak_final);
    
    fprintf(fid, 'Reconstruction (complex sum) max deviation: FIR vs delay=%.4g\n', info.recon.fir_max_abs_dev);
    fprintf(fid, 'Reconstruction (IIR explicit 3-band sum) time-domain max |low+mid+high-x| = %.4g\n', info.recon.iir_time_recon_max_err);
    fprintf(fid, 'Reconstruction (IIR residual-mid sanity) max |low+(x-low-high)+high-x| = %.4g\n', info.recon.iir_time_recon_max_err_resid);
    fprintf(fid, 'Reconstruction magnitude ripple (dB): FIR=%.3f dB, IIR=%.3f dB\n', info.recon.fir_mag_ripple_db, info.recon.iir_mag_ripple_db);
    
    c = info.fir_checks;
    v = @(x) double(x(1));
    fprintf(fid, '\nFIR sanity checks\n');
    fprintf(fid, '- Edge gains (|H(0)|, |H(pi)|):\n');
    fprintf(fid, '  Low : (%.4f, %.4f)  [expect ~1, ~0]\n', v(c.low_dc), v(c.low_nyq));
    fprintf(fid, '  Mid : (%.4f, %.4f)  [expect ~0, ~0]\n', v(c.mid_dc), v(c.mid_nyq));
    fprintf(fid, '  High: (%.4f, %.4f)  [expect ~0, ~1]\n', v(c.high_dc), v(c.high_nyq));
    fprintf(fid, '- Symmetry max|h-flip(h)| (abs/rel):\n');
    fprintf(fid, '  Low : %.3g / %.3g\n', v(c.low_sym_abs), v(c.low_sym_rel));
    fprintf(fid, '  Mid : %.3g / %.3g\n', v(c.mid_sym_abs), v(c.mid_sym_rel));
    fprintf(fid, '  High: %.3g / %.3g\n', v(c.high_sym_abs), v(c.high_sym_rel));
    fprintf(fid, '- Hsum magnitude max abs dev from 0 dB: %.4g dB\n', v(c.hsum_mag_max_abs_dev_db));
    fprintf(fid, '- Hsum group delay median: %.3f samples; max abs dev: %.3g samples\n', v(c.hsum_gd_median), v(c.hsum_gd_max_abs_dev));
    fprintf(fid, '- Best delay fit: sign=%+d, delay=%d samples (max err=%.3g, rms err=%.3g)\n', v(c.hsum_best_sign), v(c.hsum_best_delay), v(c.hsum_best_max_err), v(c.hsum_best_rms_err));
    
    fprintf(fid, 'RMS: original=%.4f, FIR+echo=%.4f, IIR+echo=%.4f\n', stats.rms_x, stats.rms_fir, stats.rms_iir);
    fprintf(fid, 'Crest: original=%.3f, FIR+echo=%.3f, IIR+echo=%.3f\n', stats.crest_x, stats.crest_fir, stats.crest_iir);
    fprintf(fid, 'Peaks after headroom: FIR=%.4f, IIR=%.4f\n', info.fir_peak_after, info.iir_peak_after);
    
    fclose(fid);
    
    % markdown report
    mdPath = fullfile(fileparts(path), 'report.md');
    fid = fopen(mdPath, 'w');
    if fid < 0, return; end
    
    fprintf(fid, '# Topic 1 Report\n\n');
    fprintf(fid, '## Setup\n');
    fprintf(fid, '- MATLAB: %s on %s\n', info.matlabVersion, info.platform);
    fprintf(fid, '- Audio: %s\n', info.audioPath);
    fprintf(fid, '- Sample rate: %.1f Hz\n\n', info.fs);
    
    fprintf(fid, '## EQ Settings\n');
    fprintf(fid, '- Crossovers: %d Hz, %d Hz\n', p.fc1, p.fc2);
    fprintf(fid, '- Gains: low=%+.1f dB, mid=%+.1f dB, high=%+.1f dB\n\n', p.gain_low_db, p.gain_mid_db, p.gain_high_db);
    
    fprintf(fid, '## Filter Info\n');
    fprintf(fid, '- FIR order: %d (delay: %d samples = %.2f ms)\n', info.fir_order, info.fir_group_delay_samples, info.fir_latency_ms);
    fprintf(fid, '- IIR: Linkwitz-Riley, orders %d/%d at fc1/fc2\n', info.iir_orders(1), info.iir_orders(2));
    fprintf(fid, '- Max pole radius: %.4f (stable)\n\n', info.iir_max_pole_radius);
    
    fprintf(fid, '## Reconstruction\n');
    fprintf(fid, '- FIR sum error: %.4g\n', info.recon.fir_max_abs_dev);
    fprintf(fid, '- IIR sum error: %.4g\n', info.recon.iir_time_recon_max_err);
    fprintf(fid, '- Magnitude ripple: FIR=%.3f dB, IIR=%.3f dB\n\n', info.recon.fir_mag_ripple_db, info.recon.iir_mag_ripple_db);
    
    fprintf(fid, '## Echo\n');
    fprintf(fid, '- Delay: %d ms, decay: %.2f, taps: %d\n\n', p.echo.delay_ms, p.echo.alpha, p.echo.numTaps);
    
    fprintf(fid, '## Output Files\n');
    fprintf(fid, '- processed_fir_eq_echo.wav\n');
    fprintf(fid, '- processed_iir_eq_echo.wav\n');
    fprintf(fid, '- Various PNG plots\n');
    
    fclose(fid);
end