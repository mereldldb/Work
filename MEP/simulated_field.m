%% Simulated field (in T) by comsol of simple Halbach magent
select30=[ 5 25; 5 25 ];% selection of 30x30 comsol matrix field
load h3comsol;                        %field from COMSOL simulation
Bz  = h3comsol(select30(1,1):select30(1,2), select30(2,1):select30(2,2));
save Bz;

