% Converts seconds -> minutes -> hours -> days

function [t, label] = secondsTime(seconds)
    mins = seconds / 60;
    hrs = mins / 60;
    dayz = hrs / 24;
    if dayz > 1
        t     = dayz;
        label = 'days';
    elseif hrs > 1
        t     = hrs;
        label = 'hrs';
    elseif mins > 1
        t     = mins;
        label = 'mins';
    else
        t     = seconds;
        label = 'secs';
    end
end