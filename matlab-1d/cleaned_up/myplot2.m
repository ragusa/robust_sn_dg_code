function myplot2(fid,phi,porder,dx,col)

figure(fid);

ncells=length(dx);
f1=reshape(phi ,porder+1,ncells);f1=f1';
x1=0.;

for i=1:ncells,
    x2=sum(dx(1:i));
    xx=[x1 x2];
    if(porder==1)
        y1=[f1(i,1),f1(i,2)];
    else
        y1=[f1(i)  ,f1(i)  ];
    end
    plot(xx,y1,col,'LineWidth',2);
    hold on; grid on
    x1=x2;
end



y1=min(min(f1)); if(y1>0),y1=0;end
y2=max(max(f1)); y2=y2*1.05;
x1=0;
x2=sum(dx);
if(length(y1>1)),y1=max(y1);end
if(length(y2>1)),y2=max(y2);end
axis([x1 x2 y1 y2])
