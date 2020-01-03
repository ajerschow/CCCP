function CCCP(cycle,fname)

%====================================================================
% CCCP v. 1.1 Complete Calculation of Coherence Pathways. 
% Copyright (C) 1998 by Alexej Jerschow.
%
% program for the simulation of coherence transfer pathways
%
% for details on the algorithm:
%
% Alexej Jerschow and Norbert Müller, J. Magn. Reson. 134 (1998) 17-29.
%====================================================================
% authors address:
%
% Alexej.Jerschow@ico.unil.ch
% Université de Lausanne
% Section de Chimie
% BCH-Dorigny, CH-1015 Lausanne
% Switzerland
% Tel: +41 21 692 3802
% Fax: +41 21 692 4035
%====================================================================
% DISCLAIMER
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% If you have not received a copy of the GNU General Public License along with
% this program write to the Free Software Foundation, Inc., 59 Temple Place -
% Suite 330, Boston, MA 02111-1307, USA.
%====================================================================
% BUG fixes:
%  - does not exit with runtime error any more when all paths
%    are canceled by too efficient phase cycle, gradients, or diffusion
%  - DispPaths can now be run autonomously
%====================================================================

%====================================================================
% TODO
% add sort/cancel step 
% after diffusion & pulsedmp
%====================================================================

clear global
global  pulses spinsys paths sigs pc receiver ...
    size_state

epsilon=1e-15;


[gamma,range,allow_paths,base,ph,pulses,recbase,grads,phrec,cutoff,epsilon,diffconst, intergrad_dels, sampledim, spin_system, J, del,abundance]=eval(cycle);



% ==================================================
% Calc phase cycle matrix

ph{length(ph)+1}=phrec;

maxlen=0;
for i=1:length(ph)
  a=ph{i};
  if(size(a,2) > maxlen)
    maxlen=size(a,2);
  end
end

pc=[];
for i=1:length(ph)
  a=ph{i};
  b=0:maxlen-1;
  b=rem(b,length(a))+1;
  pc=[pc;a(1,b)];
end
phrec=pc(length(ph),:);
pc=pc(1:size(pc,1)-1,:);

fprintf(1,'Length of phase cycle: %d\n',size(pc,2));
%==================================================

% ==================================================
% Check input

[np nc]=size(pc);

if length(base)==1
  base=base*ones(1,np);
else
  if find(sort(size(base))~=[1 np])
    fprintf(1,'base must be of size: %i %i\n', 1, np);
    return
  end
end
pc=2*pi*pc./base(ones(1,nc),:)';

if find(size(phrec)~=[1 nc])
  fprintf(1,'phrec must be of size: %i %i\n', 1, nc);
  return
end

if find(size(recbase)~=[1 1])
  fprintf(1,'recbase must be of size: 1 1\n');
  return
end
phrec=2*pi*phrec/recbase;


if find(size(pulses)~=[2 np])
  fprintf(1,'pulses must be of size: %i %i\n', 2, np);
  return
end

maxnuc=max(pulses(2,:));

if find(size(gamma)~=[1 maxnuc])
  fprintf(1,'gamma must be of size: %i %i\n', 1, maxnuc);
  return
end

sdl=length(sampledim);
gdl=size(grads);

if (gdl(1)>1 & sdl<2)
  fprintf(1,'sampledim does not contain enough elements\n');
  return
end

if ((gdl(1)==2 | gdl(1)==3) & sdl~=2)
  fprintf(1,'sampledim does not contain the right number of elements\n');
  return
end

% correct sampledim
if length(sampledim)>1
  sampledim=[sampledim(1) sampledim(2) sampledim(2)];
end

sampledim=sampledim/2;
% check allow_paths
for i=1:maxnuc
  fprintf(1,'Checking allow_paths{%i}...\n',i);
  if size(allow_paths{i},2)~=np+1
    fprintf(1,'allow_paths{%i} should have %i columns!',i,np+1);
    return
  end
end

if(gdl(2)~=np+1 & ~isempty(grads))
  fprintf(1,'grads should have %i columns!\n',np+1);
  return
end
  
if find(size(del)~=[1 np])
  fprintf(1,'del should be of size: %i %i!\n',1, np);
end


if exist('intergrad_dels') & ~isempty(intergrad_dels) & find(size(intergrad_dels)~=[1 np])
  fprintf(1,'intergrad_dels must be of size: %i %i\n', 1, np);
  return
end

switch spin_system
  case 'AX', spin_system=2; 
  case 'AMX', spin_system=3;
  case 'AMXinv', spin_system=3.01;
end
if isempty(spin_system)
  spin_system=0;
end



% ==================================================
% J evolution with pi*J*t, Hamiltonian with 2*pi*J*t*IzIz

Ixprim=.5*[0  1; 1 0];
Iyprim=.5*[0 -j; j 0];
Izprim=.5*[1  0; 0 -1];
Ipprim=Ixprim+j*Iyprim;

Ix2=kron(Ixprim,eye(2))+kron(eye(2),Ixprim);
Iy2=kron(Iyprim,eye(2))+kron(eye(2),Iyprim);
Iz2=kron(Izprim,eye(2))+kron(eye(2),Izprim);

%dw=pulses(1,:).*tan(pi/2-pulses(3,:));

switch spin_system

  case 2,
    inistate=Iz2;
    
    clear Prop;

    for i=1:size(pulses,2)
      Hp=pulses(1,i)*Iy2; %+dw(i)*Iz2;
      HJ=2*pi*J(1)*del(i)*kron(Izprim,Izprim);
      
      Prop(:,:,i)=expm(-j*Hp)*expm(-j*HJ);
      Proj1=calcord(2);
    end
    %Detect=Proj1(:,:,4);
    Detect=kron(Ipprim',eye(2));
  case 3,
    Iy3=kron(Iy2,eye(2))+kron(eye(4),Iyprim);
    Iz3=kron(Iz2,eye(2))+kron(eye(4),Izprim);

    inistate=Iz3;

    clear Prop;
    
    for i=1:size(pulses,2)
      Hp=pulses(1,i)*Iy3;%+dw(i)*Iz3;
      HJ=2*pi*del(i)*(J(1)*kron(kron(Izprim,Izprim),eye(2))+...
	  J(2)*kron(kron(Izprim,eye(2)),Izprim)+...
	  J(3)*kron(kron(eye(2),Izprim),Izprim));
      
      Prop(:,:,i)=expm(-j*Hp)*expm(-j*HJ);
      Proj1=calcord(3);
    end
    %Detect=Proj1(:,:,5);
    Detect=kron(Ipprim',eye(4));
  case 3.01
    clear Prop
    inistate=abundance(1)*kron(Iz2,eye(2))+abundance(2)*gamma(2)/gamma(1)*kron(eye(4),Izprim);
    for i=1:size(pulses,2)
      if(pulses(2,i)==1)
	Hp=pulses(1,i)*Iy2;   %+dw(i)*Iz2;
	Hp=kron(Hp,eye(2));
      else
	Hp=pulses(1,i)*Iyprim;  %+dw(i)*Izprim;
	Hp=kron(eye(4),Hp);
      end
      	
      HJ=2*pi*del(i)*(J(1)*kron(kron(Izprim,Izprim),eye(2))+J(2)*kron(kron(Izprim,eye(2)),Izprim));
      
      Prop(:,:,i)=expm(-j*Hp)*expm(-j*HJ);
    end
    clear dummy
    Proj1=calcord(2);
    for iii=1:size(Proj1,3)
      dummy(:,:,iii)=kron(Proj1(:,:,iii),ones(2));
    end
    Proj1=dummy;
    
    clear dummy;
    Proj2=calcord(1);
    for iii=1:size(Proj2,3)
      dummy(:,:,iii)=kron(ones(4),Proj2(:,:,iii));
    end
    Proj2=dummy;
    %Detect=kron(kron(Ipprim',eye(2)),eye(2));
    Detect=kron(Ipprim',eye(4));
    %Detect=kron(eye(4),Ipprim');
  otherwise
    if(spin_system~=0 & spin_system ~= 1)
      fprintf(1,'Unknown type of spin_system specified: %f\n',spin_system);  
      return
    end
end

%keyboard

% ==================================================
% Make allow_paths if necessary
% evtl. to be adapted for multinucl. case
if(0)
  np=size(pulses,2);
  if(~exist('allow_paths') | isempty(allow_paths))
    allow_paths=zeros(2*spinsys+1,np+1);
    allow_paths(spinsys+1,1)=1;
    allow_paths([spinsys spinsys+1 spinsys+2],2)=[1 1 1]';
    if np>1
      for i=2:np
	allow_paths(:,i+1)=ones(2*spinsys+1,1);
      end
    end
    allow_paths(:,np+1)=zeros(2*spinsys+1,1);
    allow_paths(spinsys+2,np+1)=1;
  end
end
% ==================================================


% ==================================================
% Make all possible paths
% Consider Weighting of sigs according to allow_paths
%   can be included later by indexing i in allow_paths(i,:) by
%   the coh. order in paths.

disp('Building possible pathways');

numnuc=max(pulses(2,:));
np=size(pulses,2);

for nuci=1:numnuc
  nucind=find(pulses(2,:)==nuci);

  fprintf(1, 'Prep Nucleus: %i: ',nuci); 
  
  paths{nuci}=zeros(1,np+1);

  for i=1:np+1
    notfirst=0;
    opaths=paths{nuci};
    spath=size(opaths,1);
    spinsys=floor((size(allow_paths{nuci},1)-1)/2);
    if((i==1) | ~isempty(find(nucind==(i-1))))
     allowi=i;
     % if i==nucind(length(nucind))+1;
     %	allowi=np+1;
     % end
	
      if(i>1 & pulses(1,i-1)==pi)
	paths{nuci}(:,i)=-paths{nuci}(:,i-1);

	retain=find(allow_paths{nuci}(spinsys+1-paths{nuci}(:,i),allowi));
	paths{nuci}=paths{nuci}(retain,:);
      else
	for ii=1:(2*spinsys+1)
	  if(allow_paths{nuci}(ii,allowi)~=0)
	    
	    if (notfirst)
	      paths{nuci}=[paths{nuci};opaths];
	    end
	    
	    paths{nuci}(size(paths{nuci},1)-spath+1:size(paths{nuci},1),i)=(-ii+spinsys+1)*ones(spath,1);
	    
	    notfirst=1;
	  end
	end
      end
    else
      fprintf(1,'use old %i %i \n', i, nuci);
      paths{nuci}(:,i)=paths{nuci}(:,i-1);
    end  
    fprintf(1, '[%i, %i] ',i,size(paths{nuci},1));
  end
  fprintf(1,'\n');

  % calc empirical damping
  % 1) due to allow_paths
  spinsys=floor((size(allow_paths{nuci},1)-1)/2);
  dummy=spinsys+1-paths{nuci}(:,[1 1+nucind]);
 
  empirsigs{nuci}=ones(size(paths{nuci},1),1);
  for tt=[1 1+nucind]
    empirsigs{nuci}=empirsigs{nuci}.*allow_paths{nuci}(spinsys+1-paths{nuci}(:,tt), tt);
  end
  
  
  % 2) due to semiideal 180 pulses
  if(spin_system==1)
    fprintf(1,'Calc. semiideal 180...\n');

    ind180=find((pulses(1,:)>3*pi/4) & (pulses(1,:)<5*pi/4) &(pulses(2,:)==nuci));
    if(ind180)
      dummy1=(paths{nuci}(:,ind180)==-paths{nuci}(:,ind180+1));
      dummy2=[];
      for cols=ind180
	dummy2=[dummy2 max(abs(paths{nuci}(:,[cols cols+1])),[],2)];
      end
      
      pulsfehler=1-pulses(ones(size(paths{nuci},1),1),ind180)/pi;
      dummy=dummy1.*(1-pulsfehler.^2/2).^dummy2+(~dummy1).*pulsfehler.^dummy2;
      empirsigs{nuci}=empirsigs{nuci}.*prod(dummy,2);
    end
  end
  
  if(0)  
    % 3) due to semiideal 90 pulses
    ind90=find((pulses(1,:)<1) & (pulses(1,:)>0) &(pulses(2,:)==nuci));
    if(ind90)
      dummy=(paths{nuci}(:,ind90) | paths{nuci}(:,ind90+1))+1-pulses(ones(size(paths{nuci},1),1),ind90);
      empirsigs{nuci}=empirsigs{nuci}.*sin(pi*prod(dummy,2)/2);
    end
  end
end


% combine paths
oallpaths=paths{1};
oempirsigs=empirsigs{1};
if numnuc>1
  for nuci=2:numnuc
    fprintf(1,'Comb %i: ', nuci);
    for i=1:size(paths{nuci},1)
      if(i==1)
	dummy=paths{nuci}(i,:);
	allpaths=oallpaths+dummy(ones(size(oallpaths,1),1),:);
	allempirsigs=oempirsigs*empirsigs{nuci}(i);
      else
	dummy=paths{nuci}(i,:);
	allpaths=[allpaths;oallpaths+dummy(ones(size(oallpaths,1),1),:)];
	allempirsigs=[allempirsigs;oempirsigs*empirsigs{nuci}(i)];
      end
      
      if(mod(i,1)==0)
	fprintf(1, '%i ',i);
      end
    end
    oallpaths=allpaths;
    oempirsigs=allempirsigs;
    fprintf(1,'\n');
  end
  clear paths
  clear empirsigs
  paths=allpaths;
  empirsigs=allempirsigs;
  clear allpaths
  clear oallpaths
  clear oempirsigs
  clear allempirsigs;

else

  allpaths=paths{1};
  allempirsigs=empirsigs{1};
  clear paths
  clear empirsigs
  paths=allpaths;
  empirsigs=allempirsigs;
  clear allpaths
  clear oallpaths
  clear oempirsigs
  clear allempirsigs
end

sizeallpaths=size(paths,1);
% ==================================================

if isempty(paths)
  disp('No pathways remaining !')
  return
end


%save paths paths;
%====================================
% Calc. Sig from PC
%====================================
% Not tested for 1-dimensional vectors
% Check if pc is constructed correctly
%====================================

disp('Calc. sig. from PC');

if (j~=sqrt(-1))
  disp('Imaginary unit redefined !');
  return
end

sigs=sum(exp(j*(-diff(paths,1,2)*pc-phrec(ones(size(paths,1),1),:))),2);


sigs=[sigs empirsigs];


% Eliminate zeroes (or close zeroes)
fprintf(1, 'Eliminating pathways, threshold: %e\n', epsilon);
retain=find(abs(prod(sigs,2))>=epsilon);

paths=paths(retain,:);
sigs=sigs(retain,:);
%normalize to PC steps
if(size(sigs,1)==0)
  fprintf(1,'\nNo pathways selected!\n');
  return
end

sigs(:,1)=sigs(:,1)/nc;
%============================
% Calculate Gradient Damping:
%============================
if (exist('grads') & ~isempty(grads))
  % calc. effgyr
  calcgrads=1;
  
  effgyr=0;
  for nuci=1:numnuc
    spinsys=floor((size(allow_paths{nuci},1)-1)/2);
    dummy=find(allow_paths{nuci}(:,1));
    effgyr=effgyr+(-dummy(1)+spinsys+1)*gamma(nuci);
  end
  
  dummy=diff(paths,1,2)*diag(gamma(pulses(2,:)));
  effgyr=cumsum([effgyr*ones(size(paths,1),1) dummy],2);
  
  % always start with 1st nucl. magnetiz. (if =0 it does not matter, anyway)
  
  
  
  disp('Calc. Gradient damping');
  %graddamping=zeros(size(paths,1),1);

  gradphase=abs(effgyr*grads'*diag(sampledim));
  
  
  if(size(gradphase,2)==3)
    gradphase=[gradphase(:,1) sqrt(gradphase(:,2).^2+gradphase(:,3).^2)];
  end
    
  
  % Z gradients
  indgz=(gradphase(:,1)<1e-60);

  indgnz1=((gradphase(:,1)<pi/2) & (gradphase(:,1)>=1e-60));
  indgnz2=(gradphase(:,1)>=pi/2);
  
  graddamping=indgz+indgnz1.*sin(indgnz1.*gradphase(:,1))./((~indgnz1)+indgnz1.*gradphase(:,1))+indgnz2./((~indgnz2)+gradphase(:,1).*indgnz2);


  % X gradients
  if(size(gradphase,2)==2)
    indgz=(gradphase(:,2)<1e-60);

    indgnz1=((gradphase(:,2)<2.44) & (gradphase(:,2)>=2.44));
    indgnz2=(gradphase(:,2)>=2.44);
  
    graddamping=[graddamping indgz+...
	indgnz1.*2.*besselj(1,indgnz1.*gradphase(:,2))./((~indgnz1)+indgnz1.*gradphase(:,2))+(2*indgnz2./((~indgnz2)+gradphase(:,2).*indgnz2)).^(1.5)/sqrt(pi)];
  end
%%%%%%%%%%%%%%%%%%%%%%%%%% (2./(kx1*xm)).^(1.5)./sqrt(pi) 
  sigsgrad=prod(graddamping,2);
  sigs=[sigs sigsgrad];

  % Eliminate by the threshold value:  cutoff
  fprintf(1, 'Eliminating pathways (after gradient), threshold: %e\n', cutoff);

  retain=find(abs(prod(sigs,2))>=cutoff);

  paths=paths(retain,:);
  effgyr=effgyr(retain,:);
  sigs=sigs(retain,:);
else
  calcgrads=0;
end

if(size(sigs,1)==0)
  fprintf(1,'\nNo pathways selected!\n');
  return
end


%==================================================
% diffusion weighting routines

if(exist('diffconst') & exist('intergrad_dels') & ~isempty(diffconst) & ...
      ~isempty(intergrad_dels) & calcgrads)

  calcdiff=1;
  
  fprintf(1, 'Calculate diffusion weighting.\n');
  for (ii=1:size(grads,1))
    if (ii==1)
      diffsigs=cumsum(effgyr*diag(grads(ii,:)),2).^2;
    else
      diffsigs=diffsigs+cumsum(effgyr*diag(grads(ii,:)),2).^2;
    end
  end
  diffsigs=diffsigs(:,1:(size(diffsigs,2)-1));
  diffsigs=exp(-diffconst*sum(diffsigs*diag(intergrad_dels),2));

  sigs=[sigs diffsigs];
  % cutoff maybe not necessary here
else
  calcdiff=0;
end

%==================================================
% Pulsedamping


if(spin_system>1)
  clear nucpaths
  if(spin_system==3.01)
    for nuci=1:max(pulses(2,:))
      notnuci=find(~(pulses(2,:)==nuci));
      spinsys=floor((size(allow_paths{nuci},1)-1)/2);
      pstart=(spinsys+1-find(allow_paths{nuci}(:,1)))*ones(size(paths,1),1);
  
      dummy=diff(paths,1,2);
      dummy(:,notnuci)=zeros(size(dummy(:,notnuci)));
      nucpaths(:,:,nuci)=[pstart cumsum(dummy,2)];
    end
  end
    
  sigpul=zeros(size(paths,1),1);
  for i=1:size(paths,1)
    state=inistate;
    
    if(spin_system==3.01)
      for ii=1:np
	state=Prop(:,:,ii)'*state*Prop(:,:,ii);
	if (pulses(2,ii)==1)
	  state=state.*Proj1(:,:,3-nucpaths(i,ii+1,1));
	else
	  state=state.*Proj2(:,:,2-nucpaths(i,ii+1,2));
	end
      end
      %sigpul(i)=abs(trace(Detect'*state));
      sigpul(i)=trace(Detect'*state);
    else
      spinsys=floor((size(allow_paths{1},1)-1)/2);
      for ii=1:np
	state=Prop(:,:,ii)'*state*Prop(:,:,ii);
	state=state.*Proj1(:,:,spinsys+1-paths(i,ii+1,1));
      end
      %sigpul(i)=abs(trace(Detect'*state));
      sigpul(i)=trace(Detect'*state);
    end
    
    fprintf(1,' [%i]',i);
  end
  fprintf(1,'\n');
  sigs(:,2)=sigs(:,2).*sigpul;
end



% Appropriate place to include convection routines

  
% ==================================================
% Sort pathways

%phases=angle(sigs)*180/pi;
%sigs=abs(sigs);

[dummy, I]=sort(abs(prod(sigs,2)));
I=flipud(I);
paths =paths(I,:);
%acqu_states=acqu_states(I,:);
sigs=sigs(I,:);
if (calcgrads)
  effgyr=effgyr(I,:);
end

if prod(sigs(1,:))==0
  disp ('All paths canceled!');
  return;
end

% save routine for debugging
% calc paths for diff nuclei
clear nucpaths
for nuci=1:max(pulses(2,:))
  notnuci=find(~(pulses(2,:)==nuci));
  spinsys=floor((size(allow_paths{nuci},1)-1)/2);
  pstart=(spinsys+1-find(allow_paths{nuci}(:,1)))*ones(size(paths,1),1);

  dummy=diff(paths,1,2);
  dummy(:,notnuci)=zeros(size(dummy(:,notnuci)));
  
  nucpaths(:,:,nuci)=[pstart cumsum(dummy,2)];
end

svar=input('Make Logfile ? > ','s');
  
saverange=1:min(500, size(paths,1));

if (calcgrads)
  effgyr=effgyr(saverange,:);
end

paths=paths(saverange,:);
sigs=sigs(saverange,:);
nucpaths=nucpaths(saverange,:,:);
if (calcgrads)
  save (svar, 'effgyr', 'paths', 'sigs', 'nucpaths', 'pulses');
else
  save (svar, 'paths', 'sigs', 'nucpaths', 'pulses');
end

% Start of Display routine

if size(paths,1)==0
  disp('NO PATHS !!!\n')
  return
end

if range(1)>size(paths,1)
  range=size(paths,1)
end

maxnuc=max(pulses(2,:));
fh=get(0, 'children');

fh=sort(get(0, 'children'));

range=range(1):min(range(length(range)),size(paths,1));
if (length(fh)>0)
  figure(fh(1));
end

clf;
DispPaths(paths(range,:),sigs(range,:),1,spinsys);
titletext=sprintf('Total CTPs');
title(titletext);

if(maxnuc>1)
  if(length(fh)>1)
    figure(fh(2));
  else
    figure;
  end
  
  clf;
  for countnucs=1:maxnuc
    subplot(maxnuc,1,countnucs);
    cla;
    DispPaths(nucpaths(range,:,countnucs),sigs(range,:),0,spinsys);
    titletext=sprintf('Nucleus %i, gamma=%e', countnucs, gamma(countnucs));
    title(titletext);
  end
end


% ==================================================
% Make Logfile
% ==================================================

%svar=input('Make Logfile ? > ','s')

if ~isempty(svar)
  np=size(pulses,2);
  svar=[svar '.log'];
  fid=fopen(svar, 'w');
  fprintf(fid, 'Name: %s;  no of steps: %i,  no of paths: %i, no. nuclei: %i\n\n', cycle, nc, sizeallpaths, maxnuc);
  fprintf(fid,'     ');

  dummy=1;
  if maxnuc>1
    dummy=0;
  end
  
  for jj=dummy:maxnuc
    for i=0:np
      fprintf(fid, '%3i', i);
    end
    fprintf(fid, '   ');
  end
  
  fprintf(fid, '   Sig       Phase   PCdmp   Puldmp   ');

  if(calcgrads)
    fprintf(fid, 'Gdmp     ');
  end
  if(calcdiff)
    fprintf(fid,'Diffdmp  ');
  end
  
  fprintf(fid,'\n');
  
  np=size(pulses,2);
  for i=1:(np+1)*(maxnuc+(maxnuc>1))+20
    fprintf(fid, '---');
  end
  fprintf(fid, '\n');

  allsigs=prod(sigs,2);
  for i=1:min(size(paths,1),500)
    fprintf(fid, '%3i) ',i);

    for ii=1:size(paths,2)
      fprintf(fid, '%3i',paths(i,ii));
    end
    fprintf(fid, '   ');
    
    if maxnuc>1
      for jj=1:maxnuc
	for ii=1:size(paths,2)
	  fprintf(fid, '%3i',nucpaths(i,ii,jj));
	end
	fprintf(fid, '   ');
      end
    end
    
    fprintf(fid, '   %7.5f %7.2f %8.5f %8.5f ', abs(allsigs(i)), ...
	  angle(allsigs(i))*180/pi,  abs(sigs(i,1)), sigs(i,2));
  
    if(calcgrads)
      fprintf(fid,'%8.5f ',sigs(i,3));
    end
    if(calcdiff)
      fprintf(fid,'%8.5f ',sigs(i,4));
    end
    fprintf(fid,'\n');
  end
end
