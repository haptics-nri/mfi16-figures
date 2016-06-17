function visforce_movie(filename, force, pose, mass, com)

    vid = VideoWriter(filename);
    vid.FrameRate = 60;
    open(vid);
    sl = visforce(force, pose, mass, com);
    
    start = get(sl, 'Min');
    stop = get(sl, 'Max');
    for i=start:50:stop
        set(sl, 'Value', i);
        feval(get(sl, 'Callback'), sl);
        drawnow;
        writeVideo(vid, getframe(gcf));
    end
    
    close(vid);

end

