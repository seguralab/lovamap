% function to associate bead size with color
% varargin:  custom color

function color = beadColor(diam, varargin)

    if ~isempty(varargin)
        color = varargin{1};
    else
        % Define colors
%         % red range
%         c40 =  '#ECDDD2';
%         c50 =  '#E8C9C0';
%         c60 =  '#E1A39B';
%         c70 =  '#DA7C75';
%         c80 =  '#D45550';
%         c90 =  '#CD2F2B';
%         c100 = '#C60806';
%         c110 = '#B90B09';
%         c120 = '#AB0D0B';
%         c130 = '#9E0F0E';
%         c140 = '#911211';
%         c150 = '#831513';
%         c160 = '#761716';
%         c170 = '#691A19';
%         c180 = '#5B1C1B';
%         c190 = '#4E1E1E';
%         c200 = '#332423';

        % gray colors range
        c40 =  [235, 236, 246] / 255;
        c50 =  [245, 246, 232] / 255;
        c60 =  [246, 234, 224] / 255;
        c70 =  [218, 234, 245] / 255;
        c80 =  [238, 231, 246] / 255;
        c90 =  [254, 236, 240] / 255;
        c100 = [231, 246, 227] / 255;
        c110 = [238, 231, 231] / 255;
        c120 = [233, 253, 248] / 255;
        c130 = [230, 234, 250] / 255;
        c140 = [254, 249, 233] / 255;
        c150 = [225, 253, 224] / 255;
        c160 = [253, 237, 248] / 255;
        c170 = [255, 233, 209] / 255;
        c180 = [219, 231, 240] / 255;
        c190 = [237, 232, 223] / 255;
        c200 = [240, 240, 240] / 255;

        % Round to nearest 10
        diam = round(diam / 10) * 10;

        if diam <= 40
            color = c40;
        elseif diam >= 200
            color = c200;
        else
            switch diam
                case 50
                    color = c50;
                case 60
                    color = c60;
                case 70
                    color = c70;    
                case 80
                    color = c80;
                case 90
                    color = c90;
                case 100
                    color = c100;
                case 110
                    color = c110;
                case 120
                    color = c120;
                case 130
                    color = c130;
                case 140
                    color = c140;
                case 150
                    color = c150;
                case 160
                    color = c160;
                case 170
                    color = c170;
                case 180
                    color = c180;
                case 190
                    color = c190;
                otherwise
                    color = c100;
            end
        end
    end
end