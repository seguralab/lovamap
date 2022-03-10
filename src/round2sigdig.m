% Option to round up or round down to the nearest significant digit
% This function was written to round data according to the dx - data would
% not be any more precise than the sig digits of dx
% You need to go through examples to see how this function rounds...I don't
% know how to properly describe the intention without examples
% varargin: 'up'
%           'down'
%           (default no input will do standard rounding)

function output = round2sigdig(number, dx, varargin)
    dec = 0;
    if dx <= 0
        error('dx must be greater than 0')
    elseif dx < 1
        while (floor(dx * 10^dec) ~= dx * 10^dec)
            dec = dec + 1;
        end
        if strcmp(varargin, 'up')
            output = round(number + (5 * 10^-dec), dec - 1);
        elseif strcmp(varargin, 'down')
            output = round(number - (5 * 10^-dec), dec - 1);
        else
            output = round(number, dec - 1);
        end
    else
        while dx / 10^dec > 1
            %mod(dx / 10^dec, 10) == 0
            dec = dec + 1;
        end
        if strcmp(varargin, 'up')
            output = round(number + (5 * 10^(dec - 1)), -dec);
        elseif strcmp(varargin, 'down')
            output = round(number - (5 * 10^(dec - 1)), -dec);
        else
            output = round(number, -dec);
        end
    end
end