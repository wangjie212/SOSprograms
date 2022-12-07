spotpath = '../../../Programs/spotless';
addpath(genpath(spotpath));

x = msspoly('x', 2);
a = [1;1];
f = [x(1)+x(2); x(1) - x(1)^2; x(2) - x(2)^2];
g = [x(1)*x(2)];
h = [x(1)^2+x(2)^2-1];
[blk, At, C, b] = sosprogram(a, f, g, h, x, 2);
sqlp(blk, At, C, b);
