% Check if array = [0]

function logical = isscalarzero(array)
    if isscalar(array) && array(1) == 0
        logical = true;
    else
        logical = false;
    end
end