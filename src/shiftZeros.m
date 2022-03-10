% Peter function
%
% Example:
% INPUT
% [1 0 5 0 0]
% OUTPUT
% [1 5 0 0 0]

function [v] = shiftZeros(v)
	len = length(v);
    if len == 0
        return;
    end

    itFast = 1;

    while v(itFast) ~= 0
        itFast = itFast + 1;
    end
    itSlow = itFast;

    while itFast <= len
        if v(itFast) ~= 0
        	v(itSlow) = v(itFast);
            v(itFast) = 0;
            itFast = itFast + 1;
            itSlow = itSlow + 1;
        else
            itFast = itFast + 1;
        end
    end
end