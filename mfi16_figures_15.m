% part 15 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
%
vv =          cell(length(materials), length(reps), length(tools));
ii =          cell(length(materials), length(reps), length(tools));
vbody =       cell(length(materials), length(reps), length(tools));
vend =        cell(length(materials), length(reps), length(tools));
vint =        cell(length(materials), length(reps), length(tools));
vbodyint =    cell(length(materials), length(reps), length(tools));
vendint =     cell(length(materials), length(reps), length(tools));
intbody =     cell(length(materials), length(reps), length(tools));
intworld =    cell(length(materials), length(reps), length(tools));
intworldsub = cell(length(materials), length(reps), length(tools));
accint =      cell(length(materials), length(reps), length(tools));
accworld =    cell(length(materials), length(reps), length(tools));
vel =         cell(length(materials), length(reps), length(tools));
