function [gamma,range,allow_paths,base,ph,pulses,recbase,grads,phrec,cutoff,epsilon,diffconst, intergrad_dels, sampledim, spin_system, J, del, abundance]=bipol()

clear ph
ph{1}= [ 0 0 1 1 2 2 3 3];
ph{2}= [ 1 2 2 3 3 0 0 1];
ph{3}= 0;
ph{4}= [ 2 2 3 3 ];
ph{5}= [ 3 3 0 0 ];
ph{6}= 2;
ph{7}= 0;
phrec = [ 0 2 0 2 2 0 2 0 ];

base=4;
recbase=4;

pulses = [ [1 2 1 1 2 1 1]*pi/2*.95;...
            1 1 1 1 1 1 1];

gamma=2.675e8;

allow_paths{1}=[0 0 1 1 1 1 1 0;...
                0 0 1 1 1 1 1 0;...
                0 1 1 1 1 1 1 0;...
                1 1 1 1 1 1 1 0;...
                0 1 1 1 1 1 1 1;...
                0 0 1 1 1 1 1 0;...
                0 0 1 1 1 1 1 0]

spin_system='AMX';
epsilon=1e-4;
cutoff=1e-4;

range=1:8;

gradunits=50e-2*1e-3;
sampledim=[1e-2 .5e-2];

grads = [.5*[0 1 -1 0 1 -1 0 0];...
        .05*[0 0  0 1 0  0 0 0];...
        .05*[0 0  0 0 0  0 1 0]]*gradunits;

del=[.002 .002 .015 .002 .002 .015 0];
J=[16 12 7];

diffconst=7e-10;
intergrad_dels=[0 .002 .002 0.05 0.002 .002 0.3];

return







