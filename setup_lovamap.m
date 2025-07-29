function setup_lovamap()
    % LOVAMAP setup script

    fprintf('--- LOVAMAP SETUP STARTING ---\n');

    %% 1. Check and install required MathWorks toolboxes
    requiredToolboxes = {
        'Image Processing Toolbox'
        'Statistics and Machine Learning Toolbox'
    };

    v = ver;
    installedToolboxes = {v.Name};

    for i = 1:length(requiredToolboxes)
        if ~ismember(requiredToolboxes{i}, installedToolboxes)
            fprintf('Installing MathWorks toolbox: %s\n', requiredToolboxes{i});
            matlab.addons.installToolbox(requiredToolboxes{i}); % Optional: Use Add-On Explorer manually if this fails
        else
            fprintf('✔ %s already installed.\n', requiredToolboxes{i});
        end
    end

    %% 2. Install third-party dependencies (if not already on path)

    externalDir = fullfile(pwd, 'external');
    if ~exist(externalDir, 'dir')
        mkdir(externalDir);
        fprintf('Created missing external/ folder. Please add dependencies.\n');
    end

    download_third_party_dependencies(externalDir);

    thirdPartyFolders = {
        'external/SLM'
        'external/circlefit3d'
        'external/inhull'
        'external/distinct_colors'
        'external/GetFullPath'
    };

    for i = 1:length(thirdPartyFolders)
        folderPath = fullfile(pwd, thirdPartyFolders{i});
        if exist(folderPath, 'dir')
            addpath(genpath(folderPath));
            fprintf('✔ Added to path: %s\n', folderPath);
        else
            warning('Dependency folder not found: %s', folderPath);
        end
    end

    %% 3. Compile MEX files
    fprintf('\n--- Checking C++ Compiler Configuration ---\n');

    if ismac
        platformNote = ' (macOS detected)';
    elseif ispc
        platformNote = ' (Windows detected)';
    elseif isunix
        platformNote = ' (Linux detected)';
    else
        platformNote = ' (Unknown OS)';
    end
    
    fprintf('Operating System: %s%s\n', computer, platformNote);

    cfg = mex.getCompilerConfigurations('C++', 'Selected');
    if isempty(cfg)
        warning('No supported C++ compiler configured.');
    
        if ismac
            fprintf(['On macOS, MATLAB requires the full Xcode installation to compile MEX files.\n' ...
                     '   Steps:\n' ...
                     '   1. Install Xcode from the App Store\n' ...
                     '   2. Run:  sudo xcodebuild -license accept\n' ...
                     '   3. In MATLAB:  mex -setup C++\n\n']);
        elseif ispc
            fprintf(['On Windows, you can install the MinGW-w64 compiler using MATLAB Add-On Explorer.\n' ...
                     '   Then run in MATLAB:  mex -setup C++\n\n']);
        elseif isunix
            fprintf('On Linux, install GCC/G++ and run:  mex -setup C++\n');
        else
            fprintf('Unknown OS — please consult MATLAB compiler setup docs.\n');
        end
    
        fprintf('Then re-run setup_lovamap to compile MEX files.\n');
        return;
    else
        fprintf('✔ Using C++ compiler: %s\n', cfg.Name);
    end
    
    try
        fprintf('Compiling MEX files...\n');
        mexDir = fullfile(pwd, 'mex-lovamap');
        cd(mexDir);
        compile_mex;
        cd('..');
        fprintf('✔ MEX files compiled successfully.\n');
    catch ME
        warning(ME.id, '%s', ME.message);
    end

    mexDir = fullfile(pwd, 'mex-lovamap');
    if exist(mexDir, 'dir')
        addpath(genpath(mexDir));
        fprintf('✔ Added source code to path: %s\n', mexDir);
    else
        warning('mex-lovamap/ folder not found. LOVAMAP may not run correctly.');
    end


    %% 4. Add /src to path
    srcDir = fullfile(pwd, 'src');
    if exist(srcDir, 'dir')
        addpath(genpath(srcDir));
        fprintf('✔ Added source code to path: %s\n', srcDir);
    else
        warning('src/ folder not found. LOVAMAP may not run correctly.');
    end


    %% 4. Save the updated path
    fprintf('--- LOVAMAP SETUP COMPLETE ---\n');
end

function download_third_party_dependencies(externalDir)
    dependencies = {
        'SLM', 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/24443/versions/15/download/zip'
        'circlefit3d', 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34792/versions/2/download/zip';
        'inhull', 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/10226/versions/3/download/zip';
        'distinct_colors', 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/29702/versions/3/download/zip';
        'GetFullPath', 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e5728cef-4a80-11e4-9553-005056977bd0/8a238c41-22be-43c2-833b-9818e66489f1/packages/zip';
    };

    for i = 1:size(dependencies, 1)
        name = dependencies{i,1};
        url = dependencies{i,2};
        targetFolder = fullfile(externalDir, name);
        if ~exist(targetFolder, 'dir')
            fprintf('⬇ Downloading and extracting %s...\n', name);
            zipFile = fullfile(externalDir, [name '.zip']);
            try
                websave(zipFile, url);
                unzip(zipFile, targetFolder);
                
                % Normalize folder contents (some zips create a nested folder)
                contents = dir(targetFolder);
                subdirs = contents([contents.isdir] & ~startsWith({contents.name}, '.'));
                if numel(subdirs) == 1
                    nested = fullfile(targetFolder, subdirs(1).name);
                    movefile(fullfile(nested, '*'), targetFolder);
                    rmdir(nested, 's');
                end

                delete(zipFile);
                fprintf('✔ %s installed in %s\n', name, targetFolder);
            catch ME
                warning('Failed to download %s: %s', name, ME.message);
            end
        else
            fprintf('✔ %s already exists. Skipping download.\n', name);
        end
    end
end