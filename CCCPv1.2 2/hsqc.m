function [gamma,range,allow_paths,base,ph,pulses,recbase,grads,phrec,...
      cutoff,epsilon,diffconst,intergrad_dels,sampledim,spin_system,J,...
      del,abundance]=hsqc()

ph{1}= [0]; ph{2}= [0]; ph{3}=0; ph{4}=1; ph{5}=[0 2]; ph{6}=[0 0 2 2];
ph{7}=[0 0 0 0 2 2 2 2]; ph{8}=0; ph{9}=ph{7}; ph{10}=0; ph{11}=0;
phrec= [0 2 0 2  2 0 2 0]; base=4; recbase=4; 

pulses =      [[1 2 2 1 1  2 2 1 1 2  2]*pi/2*.95;...
               [1 1 2 1 2  1 2 1 2 1  2]];

allow_paths{1}=[0 0 1 1 1  1 1 1 1 1  0 0;...
                0 1 1 1 1  1 1 1 1 1  0 0;...
                1 1 1 1 1  1 1 1 1 1  0 0;...
                0 1 1 1 1  1 1 1 1 1  1 1;...
                0 0 1 1 1  1 1 1 1 1  0 0];

allow_paths{2}=[0 1 1 1 1  1 1 1 1 1  1 0;...
                1 1 1 1 1  1 1 1 1 1  1 1;...
                0 1 1 1 1  1 1 1 1 1  1 0];

gamma=[2.67522e8 0.67283e8];
gradunits=50e-2*1e-3; sampledim=1e-2;
grads =    0.1*[0 0 0 0 0 0 gamma(1)/gamma(2) 0 0 0 0 1]*gradunits;

diffconst=7e-10;
intergrad_dels=[0 0 0 0 0 0       0       .0037 0 0 0];

spin_system='AMXinv'; abundance=[1 .01]; J=[16 140];
del=[0.0017 0 0.0017 0 0.015 0.017 0.002 0 0.0017 0 0.0017];

epsilon=1e-4;
cutoff=1e-4;
range=1:8;

return
