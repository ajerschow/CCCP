function [gamma,range,allow_paths,base,ph,pulses,recbase,grads,phrec,cutoff,epsilon,diffconst, intergrad_dels, sampledim, spin_system, J, del,abundance]=dalvit()

clear ph
ph{1}= [0 0 0 0  2 2 2 2]
ph{2}= 0
ph{3}= 0
ph{4}= 0
ph{5}= 0
ph{6}= 0
ph{7}= 0
ph{8}= [0 3 2 1]
phrec= [0 2 0 2  2 0 2 0]

epsilon=1e-4;
cutoff=1e-4;
gamma=2.67522e8;

J=[16 12 7];

spin_system='AMX';
range=1:8;

pulses = [ [1 1 1 2  1 2 1 2]*pi/2*.97;...
            1 1 1 1  1 1 1 1];

gradunits=1e-2*1.25e-3;
sampledim=[1e-2 .5e-2];

grads = [[0 0 0 0 0 5 -5 0 20];...
	 [0 0 -5 0 0 0  0 0 0];...
	 0*[0 0 0 2 2 0 0 0 0]]*gradunits;
 
allow_paths{1}=[0 0 1 1 1 1 1 1 0;...
                0 0 1 1 1 1 1 1 0;...
                0 1 1 1 1 1 1 1 0;...
                1 1 1 1 1 1 1 1 0;...
                0 1 1 1 1 1 1 1 1;...
                0 0 1 1 1 1 1 1 0;...
                0 0 1 1 1 1 1 1 0];
       
diffconst=7e-10;
intergrad_dels=[0 0 .050 .015 0.050 0.00135 0.00135 0.0027];

del=[[.015 .015]*1 .015 .015 .015*.5 .00135 .00135 .00135];

base=4;
recbase=4;

return









