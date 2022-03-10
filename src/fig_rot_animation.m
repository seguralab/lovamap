function fig_rot_animation(plot, filename)
%     az_init = -37.5;
%     el_init = 30;
%     az_init = 147;
%     el_init = 31;
    az_init = 122;
    el_init = 7;

    view([az_init, el_init])
    
    % frames per second
    FPS = 30;

    % These specify the starts and ends for the viewing angles for each
    % sequence. Five start/stops -> four sequences. 
%     azs = [az_init, az_init + 24.13, az_init + 109.89, az_init + 177.71, az_init];
%     els = [el_init, el_init - 31.71, el_init - 18.891, el_init - 37.386, el_init];
    azs = [az_init, az_init + 58, az_init + 112, az_init + 317, az_init];
    els = [el_init, el_init + 0, el_init + 20, el_init + 5, el_init];

    % How many seconds per sequence/transition
%     secondsPerSeq = [2.5, 3, 3, 3.7];
    secondsPerSeq = [3, 3, 3, 3];

    % Set up splines
    tArray = [0, cumsum(secondsPerSeq)];
    azCoeffs = [-2 * (azs(2) - azs(1)) / secondsPerSeq(1)^3, 3 * (azs(2) - azs(1)) / secondsPerSeq(1)^2, 0, azs(1);
                -2 * (azs(3) - azs(2)) / secondsPerSeq(2)^3, 3 * (azs(3) - azs(2)) / secondsPerSeq(2)^2, 0, azs(2);
                -2 * (azs(4) - azs(3)) / secondsPerSeq(3)^3, 3 * (azs(4) - azs(3)) / secondsPerSeq(3)^2, 0, azs(3);
                -2 * (azs(5) - azs(4)) / secondsPerSeq(4)^3, 3 * (azs(5) - azs(4)) / secondsPerSeq(4)^2, 0, azs(4)];
    elCoeffs = [-2 * (els(2) - els(1)) / secondsPerSeq(1)^3, 3 * (els(2) - els(1)) / secondsPerSeq(1)^2, 0, els(1);
                -2 * (els(3) - els(2)) / secondsPerSeq(2)^3, 3 * (els(3) - els(2)) / secondsPerSeq(2)^2, 0, els(2);
                -2 * (els(4) - els(3)) / secondsPerSeq(3)^3, 3 * (els(4) - els(3)) / secondsPerSeq(3)^2, 0, els(3);
                -2 * (els(5) - els(4)) / secondsPerSeq(4)^3, 3 * (els(5) - els(4)) / secondsPerSeq(4)^2, 0, els(4)];

    azpp = mkpp(tArray, azCoeffs);
    elpp = mkpp(tArray, elCoeffs);
    
    % Plot to see the interpolation
    % xx = linspace(tArray(1), tArray(end), 100);
    % plot(xx, ppval(azpp, xx));
    % hold on
    % plot(tArray, azs, 'o');
    % title('Path of az');
    % hold off
    % pause
    % plot(xx, ppval(elpp, xx));
    % hold on
    % plot(tArray, els, 'o');
    % title('Path of el');
    % hold off
    % pause
    
    %disp('Rendering sequence one')
    baseCount = FPS * tArray(1);
    numFrames = FPS * secondsPerSeq(1);
    xx = linspace(tArray(1), tArray(2), numFrames + 1);
    azloc = ppval(azpp, xx(1 : end - 1));
    elloc = ppval(elpp, xx(1 : end - 1));
    count = baseCount;
    for i = 1 : numFrames
        count = count + 1;
        az = azloc(i);
        el = elloc(i);
        view([az,el])
        drawnow
        f = getframe(plot);
        im{count} = frame2im(f);
    end

    %disp('Rendering sequence two')
    baseCount = FPS * tArray(2);
    numFrames = FPS * secondsPerSeq(2);
    xx = linspace(tArray(2), tArray(3), numFrames + 1);
    azloc = ppval(azpp, xx(1 : end - 1));
    elloc = ppval(elpp, xx(1 : end - 1));
    count = baseCount;
    for i = 1 : numFrames
        count = count + 1;
        az = azloc(i);
        el = elloc(i);
        view([az,el])
        drawnow
        f = getframe(plot);
        im{count} = frame2im(f);
    end    

    %disp('Rendering sequence three')
    baseCount = FPS * tArray(3);
    numFrames = FPS * secondsPerSeq(3);
    xx = linspace(tArray(3), tArray(4), numFrames + 1);
    azloc = ppval(azpp, xx(1 : end - 1));
    elloc = ppval(elpp, xx(1 : end - 1));
    count = baseCount;
    for i = 1 : numFrames
        count = count + 1;
        az = azloc(i);
        el = elloc(i);
        view([az,el])
        drawnow
        f = getframe(plot);
        im{count} = frame2im(f);
    end  

    %disp('Rendering sequence four')
    baseCount = FPS * tArray(4);
    numFrames = FPS * secondsPerSeq(4);
    xx = linspace(tArray(4), tArray(5), numFrames + 1);
    azloc = ppval(azpp, xx);
    elloc = ppval(elpp, xx);
    count = baseCount;
    for i = 1 : numFrames + 1
        count = count + 1;
        az = azloc(i);
        el = elloc(i);
        view([az,el])
        drawnow
        f = getframe(plot);
        im{count} = frame2im(f);
    end

    file_name = filename;
    for idx = 1:count
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,file_name,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,file_name,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end