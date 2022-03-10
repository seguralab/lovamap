function [] = writeTimesToFile(lovamapData, timeLog, outputFile)

    fileHandle = fopen(outputFile, 'wt+');

    % Extract all the timing names into a cell array
    timeDataAsCell = struct2cell(timeLog(:));
    timeDataNames = timeDataAsCell(1, :);

    % Find the name of longest length
    outputNames = {'Timestamp'; 'LINEBREAK'};
    outputValues = {lovamapData.Timestamp; 'LINEBREAK'};

    outputNames = [outputNames; fieldnames(lovamapData.InputArgs)];

    fnames = fieldnames(lovamapData.InputArgs);
    for i = 1 : length(fnames)
        outputValues = [outputValues; num2str(lovamapData.InputArgs.(fnames{i}))];
    end

    outputNames = [outputNames; 'LINEBREAK'];
    outputValues = [outputValues; 'LINEBREAK'];

    outputNames = [outputNames; 'Number of voxels'];
    outputValues = [outputValues; num2str(lovamapData.numVoxels)];

    outputNames = [outputNames; 'Number of beads'];
    outputValues = [outputValues; num2str(length(lovamapData.beads))];

    outputNames = [outputNames; 'LINEBREAK'];
    outputValues = [outputValues; 'LINEBREAK'];

    outputNames = [outputNames; timeDataNames(:)];

    for i = 1 : length(timeLog)
        outputValues = [outputValues; compose('%.5f %s', timeLog(i).Time, timeLog(i).Units)];
    end

    outputNames = pad(outputNames, 'left');

    for i = 1 : length(outputNames)
        if contains(outputNames{i}, 'LINEBREAK')
            fprintf(fileHandle, '\n');
        else
            fprintf(fileHandle, '%s: %s\n', outputNames{i}, outputValues{i});
        end
    end

    fclose(fileHandle);
end
