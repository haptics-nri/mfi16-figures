% part 51 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
    transformed_grav(i,2:4) = freecalib.R * freecalib.ideal_grav(i,2:4)' + freecalib.t;
