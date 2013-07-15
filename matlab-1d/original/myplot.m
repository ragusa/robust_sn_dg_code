function myplot(fid,phin,phi,porder,dx)

figure(fid);

ncells=length(dx);
f1=reshape(phin,porder+1,ncells);f1=f1';
f2=reshape(phi ,porder+1,ncells);f2=f2';
x1=0.;

for i=1:ncells,
    x2=sum(dx(1:i));
    xx=[x1 x2];
    if(porder==1)
        y1=[f1(i,1),f1(i,2)];
        y2=[f2(i,1),f2(i,2)];
    else
        y1=[f1(i)  ,f1(i)  ];
        y2=[f2(i)  ,f2(i)  ];
    end
    plot(xx,y1,'r+-',xx,y2,'bo-','LineWidth',2);
    hold on; grid on
    x1=x2;
end

legend('Reduced Upwind','Standard Upwind','Location','Best');
xlabel('position','FontSize',12);
ylabel('Scalar flux','FontSize',12);

y1=min(min(f1,f2)); if(y1>0),y1=0;,end
y2=max(max(f1,f2)); y2=y2*1.05;
x1=0;
x2=sum(dx);
if(length(y1>1)),y1=max(y1);end
if(length(y2>1)),y2=max(y2);end
axis([x1 x2 y1 y2])
