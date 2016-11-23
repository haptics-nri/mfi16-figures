% part 9 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% load in some manually picked values

% NB to find these values:
% 1. plot(d.v(:,1)-d.v(1,1), d.v(:,2:4), d.int(:,1)-d.v(1,1), d.int(:,2:4))
% 2. d.off = -mean([force spike maxima] - [position spike minima]);
% 3. run the icra17_process steps
% 4. plot(d.iws(:,2:4))
% 5. d.ss = [index just after second spike, index just before third spike];
% 6. plot(d.biws(:,2:4))
% 7. d.bss = [index just after second spike, index just before third spike];

d = data14('abs');
    d.off = 16.0213;
    d.ss =  [28490 125800];
    d.bss = [15400 113000];
data14('abs') = d;
d = data14('glitter');
    d.off = 16.0197;
    d.ss =  [32500 121400];
    d.bss = [20130 108500];
data14('glitter') = d;
d = data14('silk');
    d.off = 16.0192;
    d.ss =  [31560 124600];
    d.bss = [18710 111500];
data14('silk') = d;
d = data14('vinyl');
    d.off = 16.0197;
    d.ss =  [27220 131200];
    d.bss = [14260 118600];
data14('vinyl') = d;
d = data14('wood');
    d.off = 16.0158;
    d.ss =  [27520 128800];
    d.bss = [14210 116000];
data14('wood') = d;

d = data38('abs');
    d.off = 15.9987;
    d.ss =  [24647 140441];
    d.bss = [11935 127889];
data38('abs') = d;
d = data38('glitter');
    d.off = 16.0027;
    d.ss =  [30130 124500];
    d.bss = [17800 111800];
data38('glitter') = d;
d = data38('silk');
    d.off = 16.0005;
    d.ss =  [28960 121400];
    d.bss = [17220 108900];
data38('silk') = d;
d = data38('vinyl');
    d.off = 16.0048;
    d.ss =  [26920 130700];
    d.bss = [14470 118200];
data38('vinyl') = d;
d = data38('wood');
    d.off = 16.0030;
    d.ss =  [29670 124600];
    d.bss = [16840 111800];
data38('wood') = d;

