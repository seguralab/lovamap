% mex -setup C++
% 
% mex CXXOPTIMFLAGS='-O3 -fwrapv -DNDEBUG' ...
%     CXXDEBUGFLAGS='' ...
%     -silent -I./include/ ./sssr_mex.cpp


% Get current platform info
isMac = ismac;
isWindows = ispc;
isLinux = isunix && ~ismac;

% Common flags
commonFlags = {
    'CXXOPTIMFLAGS=-O3 -fwrapv -DNDEBUG', ...
    'CXXDEBUGFLAGS=', ...
    '-I./include', ...
    './sssr_mex.cpp'
};

if isMac
    % macOS-specific rpath fix
    engineLibPath = fullfile(matlabroot, 'extern', 'bin', 'maci64');
    mex('-silent', ...
        commonFlags{:}, ...
        ['-L' engineLibPath], ...
        '-lMatlabEngine', ...
        '-lMatlabDataArray', ...
        ['LDFLAGS=$LDFLAGS -Wl,-rpath,' engineLibPath]);
elseif isWindows
    % Windows (assumes MinGW is properly configured)
    mex('-silent', ...
        commonFlags{:});
elseif isLinux
    % Linux: usually just works with standard flags
    mex('-silent', ...
        commonFlags{:});
else
    error('Unsupported OS for MEX compilation.');
end