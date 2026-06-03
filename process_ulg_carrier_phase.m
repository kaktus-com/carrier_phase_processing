function carrierPhaseTable = process_ulg_carrier_phase(ulgFile, varargin)
%PROCESS_ULG_CARRIER_PHASE Extract and plot carrier phase from a PX4 ULG file.
%
% Usage:
%   carrierPhaseTable = process_ulg_carrier_phase("log_135_2026-5-22-11-46-56.ulg");
%   carrierPhaseTable = process_ulg_carrier_phase("log_135_2026-5-22-11-46-56.ulg", ...
%       "Plot", true, "RelativePhase", true, "Meters", false, "HighPassHz", 5);

    arguments
        ulgFile {mustBeTextScalar} = ""
    end

    arguments (Repeating)
        varargin
    end

    opts = parseOptions(varargin{:});

    if strlength(string(ulgFile)) == 0
        [fileName, pathName] = uigetfile( ...
            {"*.ulg", "PX4 ULog files"; "*.*", "All files"}, ...
            "Select ULG file");

        if isequal(fileName, 0)
            error("No ULG file selected.");
        end

        ulgFile = fullfile(pathName, fileName);
    end

    ulgFile = string(ulgFile);
    if ~isfile(ulgFile)
        error("ULG file not found: %s", ulgFile);
    end

    scriptDir = string(fileparts(mfilename("fullpath")));
    helperScript = fullfile(scriptDir, "extract_ulg_carrier_phase.py");
    if ~isfile(helperScript)
        error("Missing helper script: %s", helperScript);
    end

    outputCsv = opts.OutputCsv;
    if strlength(outputCsv) == 0
        [inputDir, inputName] = fileparts(ulgFile);
        outputCsv = fullfile(inputDir, inputName + "_carrier_phase.csv");
    end

    pythonExe = findPython(scriptDir);
    command = sprintf('%s %s %s --output %s', ...
        shellQuote(pythonExe), shellQuote(helperScript), shellQuote(ulgFile), shellQuote(outputCsv));

    if opts.KeepZero
        command = command + " --keep-zero";
    end

    fprintf("Extracting carrier phase from %s\n", ulgFile);
    [status, output] = system(command);
    if status ~= 0
        error("Carrier phase extraction failed:\n%s", output);
    end

    carrierPhaseTable = readtable(outputCsv, "TextType", "string");
    carrierPhaseTable = addHighPassCarrierPhase(carrierPhaseTable, opts.HighPassHz);
    writetable(carrierPhaseTable, outputCsv);
    fprintf("Saved %d carrier-phase samples to %s\n", height(carrierPhaseTable), outputCsv);

    if opts.Plot
        plotCarrierPhase(carrierPhaseTable, ulgFile, opts.RelativePhase, opts.Meters, opts.HighPassHz);
    end
end

function opts = parseOptions(varargin)
    parser = inputParser;
    addParameter(parser, "OutputCsv", "", @(value) ischar(value) || isstring(value));
    addParameter(parser, "Plot", true, @(value) islogical(value) || isnumeric(value));
    addParameter(parser, "RelativePhase", true, @(value) islogical(value) || isnumeric(value));
    addParameter(parser, "Meters", false, @(value) islogical(value) || isnumeric(value));
    addParameter(parser, "KeepZero", false, @(value) islogical(value) || isnumeric(value));
    addParameter(parser, "HighPassHz", 5, @(value) isnumeric(value) && isscalar(value) && value > 0);
    parse(parser, varargin{:});

    opts = parser.Results;
    opts.OutputCsv = string(opts.OutputCsv);
    opts.Plot = logical(opts.Plot);
    opts.RelativePhase = logical(opts.RelativePhase);
    opts.Meters = logical(opts.Meters);
    opts.KeepZero = logical(opts.KeepZero);
    opts.HighPassHz = double(opts.HighPassHz);
end

function carrierPhaseTable = addHighPassCarrierPhase(carrierPhaseTable, cutoffHz)
    carrierPhaseTable.carrier_phase_cycles_highpass = NaN(height(carrierPhaseTable), 1);
    carrierPhaseTable.carrier_phase_meters_highpass = NaN(height(carrierPhaseTable), 1);

    satellites = unique(carrierPhaseTable.satellite, "stable");
    skippedSatellites = strings(0, 1);
    sampleRatesHz = zeros(0, 1);

    for index = 1:numel(satellites)
        satellite = satellites(index);
        satelliteRows = find(carrierPhaseTable.satellite == satellite);
        [timeS, order] = sort(carrierPhaseTable.time_s(satelliteRows));
        sortedRows = satelliteRows(order);
        phaseCycles = carrierPhaseTable.carrier_phase_cycles(sortedRows);

        [filteredCycles, sampleRateHz, didFilter] = highPassOneSatellite(timeS, phaseCycles, cutoffHz);
        if didFilter
            carrierPhaseTable.carrier_phase_cycles_highpass(sortedRows) = filteredCycles;
            carrierPhaseTable.carrier_phase_meters_highpass(sortedRows) = ...
                filteredCycles * 0.190293672798;
        else
            skippedSatellites(end + 1, 1) = satellite; %#ok<AGROW>
            sampleRatesHz(end + 1, 1) = sampleRateHz; %#ok<AGROW>
        end
    end

    if ~isempty(skippedSatellites)
        maxRateHz = max(sampleRatesHz, [], "omitnan");
        warning("High-pass filter was skipped for %d satellites. A %.3g Hz cutoff needs sample rate > %.3g Hz; max detected rate was %.3g Hz.", ...
            numel(skippedSatellites), cutoffHz, 2 * cutoffHz, maxRateHz);
    end
end

function [filteredPhase, sampleRateHz, didFilter] = highPassOneSatellite(timeS, phase, cutoffHz)
    filteredPhase = NaN(size(phase));
    didFilter = false;

    valid = isfinite(timeS) & isfinite(phase);
    if nnz(valid) < 4
        sampleRateHz = NaN;
        return;
    end

    validTime = timeS(valid);
    validPhase = phase(valid);
    dt = median(diff(unique(validTime)));
    if ~isfinite(dt) || dt <= 0
        sampleRateHz = NaN;
        return;
    end

    sampleRateHz = 1 / dt;
    if sampleRateHz <= 2 * cutoffHz
        return;
    end

    uniformTime = (validTime(1):dt:validTime(end)).';
    if numel(uniformTime) < 4
        return;
    end

    uniformPhase = interp1(validTime, validPhase, uniformTime, "linear", "extrap");
    uniformPhase = uniformPhase - mean(uniformPhase, "omitnan");

    sampleCount = numel(uniformPhase);
    spectrum = fft(uniformPhase);
    frequencies = (0:sampleCount - 1).' * sampleRateHz / sampleCount;
    highPassMask = frequencies >= cutoffHz & frequencies <= sampleRateHz - cutoffHz;
    spectrum(~highPassMask) = 0;
    filteredUniform = real(ifft(spectrum));

    filteredPhase(valid) = interp1(uniformTime, filteredUniform, validTime, "linear", NaN);
    didFilter = true;
end

function plotCarrierPhase(carrierPhaseTable, ulgFile, relativePhase, plotMeters, cutoffHz)
    if plotMeters
        phaseColumn = "carrier_phase_meters";
        filteredPhaseColumn = "carrier_phase_meters_highpass";
        unitsLabel = "meters";
    else
        phaseColumn = "carrier_phase_cycles";
        filteredPhaseColumn = "carrier_phase_cycles_highpass";
        unitsLabel = "cycles";
    end

    satellites = unique(carrierPhaseTable.satellite, "stable");

    figure("Name", "ULG Carrier Phase High-Pass", "Color", "w");
    tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");

    rawAxis = nexttile;
    plotPhasePanel(rawAxis, carrierPhaseTable, satellites, phaseColumn, relativePhase);
    title(rawAxis, "Before high-pass", "Interpreter", "none");

    filteredAxis = nexttile;
    plotPhasePanel(filteredAxis, carrierPhaseTable, satellites, filteredPhaseColumn, false);
    title(filteredAxis, sprintf("After %.3g Hz high-pass", cutoffHz), "Interpreter", "none");

    xlabel(rawAxis, "Time since first sensor_gps_raw sample (s)");
    xlabel(filteredAxis, "Time since first sensor_gps_raw sample (s)");

    if relativePhase
        ylabel(rawAxis, "Carrier phase relative to first sample (" + unitsLabel + ")");
    else
        ylabel(rawAxis, "Carrier phase (" + unitsLabel + ")");
    end
    ylabel(filteredAxis, "High-pass carrier phase (" + unitsLabel + ")");

    if all(~isfinite(carrierPhaseTable.(filteredPhaseColumn)))
        text(filteredAxis, 0.5, 0.5, ...
            sprintf("No filtered data\\n%.3g Hz cutoff needs sample rate > %.3g Hz", cutoffHz, 2 * cutoffHz), ...
            "Units", "normalized", "HorizontalAlignment", "center");
    end

    sgtitle("Carrier Phase from " + string(getFileName(ulgFile)), "Interpreter", "none");
end

function plotPhasePanel(axisHandle, carrierPhaseTable, satellites, phaseColumn, relativePhase)
    hold(axisHandle, "on");
    grid(axisHandle, "on");

    plottedAny = false;
    for index = 1:numel(satellites)
        satellite = satellites(index);
        rows = carrierPhaseTable.satellite == satellite;
        timeS = carrierPhaseTable.time_s(rows);
        phase = carrierPhaseTable.(phaseColumn)(rows);

        if relativePhase
            firstValid = find(isfinite(phase), 1, "first");
            if ~isempty(firstValid)
                phase = phase - phase(firstValid);
            end
        end

        if any(isfinite(phase))
            plot(axisHandle, timeS, phase, ".-", "DisplayName", char(satellite));
            plottedAny = true;
        end
    end

    if plottedAny
        legend(axisHandle, "Location", "eastoutside");
    end

    hold(axisHandle, "off");
end

function pythonExe = findPython(scriptDir)
    candidates = [
        fullfile(scriptDir, "venv", "bin", "python")
        fullfile(scriptDir, ".venv", "bin", "python")
        "python3"
        "python"
    ];

    for index = 1:numel(candidates)
        candidate = candidates(index);
        if isfile(candidate)
            pythonExe = candidate;
            return;
        end

        [status, ~] = system("command -v " + shellQuote(candidate));
        if status == 0
            pythonExe = candidate;
            return;
        end
    end

    error("No Python executable found. Install pyulog or use the local venv.");
end

function quoted = shellQuote(value)
    value = char(string(value));
    quoted = "'" + string(strrep(value, "'", "'\''")) + "'";
end

function fileName = getFileName(pathName)
    [~, name, ext] = fileparts(pathName);
    fileName = string(name) + string(ext);
end
