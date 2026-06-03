function obsTable = processing_t(rinexFile)
%PROCESSING_T Open a RINEX observation file and plot GPS L1C carrier phase.
%
% This function is self-contained so it does not require MATLAB toolboxes.
% It reads receiver.obs or receiver.txt from this folder by default.
%
% Usage:
%   obsTable = processing_t("receiver.obs");
%   obsTable = processing_t("my_rinex_file.txt");
%   obsTable = processing_t();

    arguments
        rinexFile {mustBeTextScalar} = ""
    end

    close all; clc;

    scriptDir = fileparts(mfilename("fullpath"));
    if strlength(string(rinexFile)) > 0
        obsFile = string(rinexFile);
    else
        obsFile = fullfile(scriptDir, "receiver.obs");
        if ~isfile(obsFile) && isfile(fullfile(scriptDir, "receiver.txt"))
            obsFile = fullfile(scriptDir, "receiver.txt");
        end
    end
    carrierObsType = "L1C";

    plotRelativePhase = true;   % true: subtract first valid phase per satellite
    plotInMeters = false;       % true: convert L1 cycles to meters
    l1WavelengthMeters = 0.190293672798;

    if ~isfile(obsFile)
        [fileName, pathName] = uigetfile( ...
            {"*.obs;*.rnx;*.txt;*.??o", "RINEX observation files"; "*.*", "All files"}, ...
            "Select RINEX observation file");

        if isequal(fileName, 0)
            error("No observation file selected.");
        end

        obsFile = fullfile(pathName, fileName);
    end

    fprintf("Reading %s\n", obsFile);
    obsTable = readRinexCarrierPhase(obsFile, carrierObsType);

    if isempty(obsTable)
        error("No %s carrier-phase observations were found in %s.", carrierObsType, obsFile);
    end

    fprintf("Loaded %d %s carrier-phase samples from %d satellites.\n", ...
        height(obsTable), carrierObsType, numel(unique(obsTable.Satellite)));

    figure("Name", carrierObsType + " Carrier Phase", "Color", "w");
    hold on;
    grid on;

    satellites = unique(obsTable.Satellite, "stable");
    for k = 1:numel(satellites)
        sat = satellites(k);
        satRows = obsTable.Satellite == sat;
        phase = obsTable.CarrierPhaseCycles(satRows);

        if plotRelativePhase
            firstValid = find(isfinite(phase), 1, "first");
            if ~isempty(firstValid)
                phase = phase - phase(firstValid);
            end
        end

        if plotInMeters
            phase = phase * l1WavelengthMeters;
        end

        plot(obsTable.Time(satRows), phase, ".-", "DisplayName", char(sat));
    end

    xlabel("GPS time");
    if plotInMeters
        unitsLabel = "meters";
    else
        unitsLabel = "cycles";
    end

    if plotRelativePhase
        ylabel(carrierObsType + " carrier phase relative to first sample (" + unitsLabel + ")");
    else
        ylabel(carrierObsType + " carrier phase (" + unitsLabel + ")");
    end

    title(carrierObsType + " Carrier Phase from " + string(getFileName(obsFile)), "Interpreter", "none");
    legend("Location", "eastoutside");
    hold off;
end

function obsTable = readRinexCarrierPhase(obsFile, carrierObsType)
    fid = fopen(obsFile, "r");
    if fid < 0
        error("Could not open %s.", obsFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    obsTypesBySystem = containers.Map("KeyType", "char", "ValueType", "any");

    while true
        line = fgetl(fid);
        if ~ischar(line)
            error("Reached end of file before END OF HEADER.");
        end

        label = rinexLabel(line);
        if strcmp(label, "SYS / # / OBS TYPES")
            systemId = line(1);
            if isKey(obsTypesBySystem, systemId)
                obsTypes = obsTypesBySystem(systemId);
            else
                obsTypes = {};
            end

            obsTypeText = extractColumns(line, 8, 60);
            newTypes = regexp(strtrim(obsTypeText), "\s+", "split");
            newTypes = newTypes(~cellfun("isempty", newTypes));
            obsTypesBySystem(systemId) = [obsTypes, newTypes];
        elseif strcmp(label, "END OF HEADER")
            break;
        end
    end

    times = NaT(0, 1, "TimeZone", "UTC");
    satellites = strings(0, 1);
    carrierPhaseCycles = zeros(0, 1);

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        if isempty(strtrim(line)) || line(1) ~= ">"
            continue;
        end

        epochParts = regexp(strtrim(line(2:end)), "\s+", "split");
        if numel(epochParts) < 8
            continue;
        end

        epochTime = makeEpochTime(epochParts);
        epochFlag = str2double(epochParts{7});
        numSatellites = str2double(epochParts{8});

        if isnan(numSatellites)
            continue;
        end

        for satIndex = 1:numSatellites
            obsLine = fgetl(fid);
            if ~ischar(obsLine) || strlength(string(obsLine)) < 3
                continue;
            end

            satelliteId = string(strtrim(extractColumns(obsLine, 1, 3)));
            if strlength(satelliteId) < 2
                continue;
            end

            systemId = char(extractBetween(satelliteId, 1, 1));
            if ~isKey(obsTypesBySystem, systemId)
                continue;
            end

            obsTypes = string(obsTypesBySystem(systemId));
            phaseIndex = find(obsTypes == carrierObsType, 1);
            if isempty(phaseIndex)
                phaseIndex = find(startsWith(obsTypes, "L"), 1);
            end

            if isempty(phaseIndex) || ~(epochFlag == 0 || epochFlag == 1)
                continue;
            end

            phaseCycles = readObservationValue(obsLine, phaseIndex);
            if ~isfinite(phaseCycles)
                continue;
            end

            times(end + 1, 1) = epochTime; %#ok<AGROW>
            satellites(end + 1, 1) = satelliteId; %#ok<AGROW>
            carrierPhaseCycles(end + 1, 1) = phaseCycles; %#ok<AGROW>
        end
    end

    obsTable = table( ...
        times, categorical(satellites), carrierPhaseCycles, ...
        carrierPhaseCycles * 0.190293672798, ...
        'VariableNames', {'Time', 'Satellite', 'CarrierPhaseCycles', 'CarrierPhaseMeters'});
end

function value = readObservationValue(obsLine, observationIndex)
    fieldStart = 4 + (observationIndex - 1) * 16;
    fieldEnd = fieldStart + 13;

    if strlength(string(obsLine)) < fieldStart
        value = NaN;
        return;
    end

    valueText = strtrim(extractColumns(obsLine, fieldStart, fieldEnd));
    if valueText == ""
        value = NaN;
    else
        value = str2double(valueText);
    end
end

function epochTime = makeEpochTime(epochParts)
    year = str2double(epochParts{1});
    month = str2double(epochParts{2});
    day = str2double(epochParts{3});
    hour = str2double(epochParts{4});
    minute = str2double(epochParts{5});
    second = str2double(epochParts{6});

    epochTime = datetime(year, month, day, hour, minute, floor(second), ...
        "TimeZone", "UTC") + seconds(second - floor(second));
end

function label = rinexLabel(line)
    if strlength(string(line)) < 61
        label = "";
    else
        label = char(strtrim(extractColumns(line, 61, strlength(string(line)))));
    end
end

function text = extractColumns(line, firstColumn, lastColumn)
    line = char(line);
    if numel(line) < firstColumn
        text = "";
        return;
    end

    lastColumn = min(lastColumn, numel(line));
    text = string(line(firstColumn:lastColumn));
end

function fileName = getFileName(pathName)
    [~, name, ext] = fileparts(pathName);
    fileName = string(name) + string(ext);
end
