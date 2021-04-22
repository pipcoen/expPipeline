
load X
addpath('Isomap');
[Xstrain,XGraph] = Run3dIsomap(X);
save Xstrain Xstrain
save XGraph XGraph