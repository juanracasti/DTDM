
Deltar = 7;
v = [.01 .02 .05 .1 .2 .5 1 2 5 10 50 100];

%Load H and theta and interpolate

filename = ['H_Dr_' num2str(Deltar) '.mat'];
load(filename);

g_grid = 0:.002:1;
g_itp  = 0:.0001:1;

H_itp = interp2(g_grid,theta,H,g_itp,theta,'spline');

for i = 1:length(v)
    
    new_runchrono(H_itp,theta,Deltar,v(i),1)
    
end

for i = 4:length(v)
    
    new_runchrono(H_itp,theta,Deltar,1,v(i))
    
end