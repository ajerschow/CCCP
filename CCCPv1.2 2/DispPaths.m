function DispPaths(paths,sigs,drawleg,spinsys)

%global spinsys

maxorder=max(3,spinsys+1); 

np=size(paths,2);

index=[1:size(paths,2)];
index=[index;index];
index2=reshape(index,1,2*size(paths,2));

paths=paths(:,index2);
elev=.05*([1:size(paths,1)]'-size(paths,1)/2);
elev=elev(:,ones(1,size(paths,2)));
paths=paths+elev;

index(1,:)=index(1,:)-.1;
index(2,:)=index(2,:)+.1;
indexsize=size(index,2);
index2=reshape(index,1,2*indexsize);
index2=[0 index2];
index2=index2(1:2*indexsize);

hh=plot(index2(ones(size(paths,1),1),:)',paths');
h=gca;
set(h,'XTick',[0:np]);
set(h,'YTick',[-maxorder:maxorder]);
set(h,'YLim',[-maxorder maxorder]);
set(h,'XLim',[0 np]);
set(h,'XGrid','on');
set(h,'YGrid','on');

allsigs=prod(sigs,2);
for i=1:length(hh)
  set(hh(i),'LineWidth',10*0.5*abs(allsigs(i))/abs(allsigs(1))+1e-15);
end

if(drawleg)
  for i=1:size(paths,1)
    out=sprintf('%5.4e %4.3f ph %3i', abs(allsigs(i)), abs(allsigs(i))/abs(allsigs(1)), round(angle(allsigs(i))*180/pi));
    if (i==1)
      leg=out;
    else
      leg=str2mat(leg,out);
    end
  end  

  leg=legend(leg,-1);
end

%----------------------------







