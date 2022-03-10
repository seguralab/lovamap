% Option to round up or round down to the nearest significant digit
% varargin: 'up'
%           'down'
%           (default no input will do standard rounding)

function output = round2sigdig_exact(number, dx, varargin)
    dec = 0;
    inte = 0;
    if number == dx
        output = number;
    elseif dx <= 0
        error('dx must be greater than 0')
    elseif dx < 1
        while (floor(dx * 10^dec) ~= dx * 10^dec)
            dec = dec + 1;
        end
        if strcmp(varargin, 'up')
            output = round(number + (5 * 10^-(dec + 1)), dec);
        elseif strcmp(varargin, 'down')
            output = round(number - (5 * 10^-(dec + 1)), dec);
        else
            output = round(number, dec);
        end
    else
        while (floor(dx * 10^dec) ~= dx * 10^dec)
            dec = dec + 1;
        end
        while dx / 10^inte > 1
            inte = inte + 1;
        end
        if dec > 0
            if strcmp(varargin, 'up')
                output = round(number + (5 * 10^-(dec + 1)), dec);
            elseif strcmp(varargin, 'down')
                output = round(number - (5 * 10^-(dec + 1)), dec);
            else
                output = round(number, dec);
            end
        else
            if strcmp(varargin, 'up')
                output = round(number + (5 * 10^(inte - 2)), -(inte - 1));
            elseif strcmp(varargin, 'down')
                output = round(number - (5 * 10^(inte - 2)), -(inte - 1));
            else
                output = round(number, -(inte - 1));
            end
        end
    end
end