clear all

% set bounds (in degrees and minutes)
% [ minlat� minlat' maxlat� maxlat' minlon� minlon' maxlon� maxlon']
% l = [40 00 52 00  -5 00   8 00]; % France
% l = [43 48 44 05   4 30   5 00]; % Cadarache, France
l = [38 00 39 00  20 00  21 00]; % Kefalonia, Greece
% l = [37 10 37 40 138 15 138 55]; % Kashiwazaki, Japan
% l = [18 00 21 50 -158 50 -154 00]; % Mauna Loa, Hawai

% transform coordinates
latcrop = [l(1)+l(2)/60 l(3)+l(4)/60];
loncrop = [l(5)+l(6)/60 l(7)+l(8)/60];

% equivalent sizes
[~,xlat] = lonlat2m([1 1]*mean(loncrop),latcrop);
xlon = lonlat2m(loncrop,[1 1]*mean(latcrop));
s = abs([diff(xlon) diff(xlat)])/1000

