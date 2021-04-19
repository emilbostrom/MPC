function rita(X,U)
global xstar
if isempty(X)
    msgbox('Simulate first!','Help','help')
    return
  end
  c =[1 0 0;0 0.65 0;0 0 1;0 0 0;0.8 0.8 0;0.75 0 0.75;0.5 0.5 0.8];
  n=1;
  figure(1)
  while (isempty(findobj(gcf,'color',c(n,:)))==0)&(n<length(c))
    n=n+1;
  end
  if length(X.signals)==2
    t=X.time;
    z1=X.signals(1).values;
    z2=X.signals(2).values;   
  else
    t=X.time;
    z1=X.signals.values(:,1);
    z2=X.signals.values(:,8);
  end
  subplot(211)
  plot(t,z1,'color',c(n,:)), hold on
  if length(X.signals)==2
     h1=get(gca,'Children');
     set(h1(1),'Linestyle','-.')
  end
  plot(t,xstar(1)*ones(size(t)),'k:','linewidth',3);
  xlabel('Time [min]','fontsize',14)
  ylabel('z1','fontsize',14)
  title('Measured signals','fontsize',14)
  axis tight
  yl=ylim;
  axis([xlim yl+[-1 1]*0.05*diff(yl)])

  subplot(212)
  plot(t,z2,'color',c(n,:)), hold on
  if length(X.signals)==2
     h1=get(gca,'Children');
     set(h1(1),'Linestyle','-.')
  end
  plot(t,xstar(8)*ones(size(t)),'k:','linewidth',3);
  xlabel('Time [min]','fontsize',14)
  ylabel('z2','fontsize',14)
  axis tight
  yl=ylim;
  axis([xlim yl+[-1 1]*0.05*diff(yl)])

  
  figure(2)
  t=U.time;
  u1=U.signals(1).values;
  u2=U.signals(2).values;
  subplot(211)
  plot(t,u1,'color',c(n,:)), hold on
  xlabel('Time [min]','fontsize',14)
  ylabel('V','fontsize',14)
  title('Control signals','fontsize',14)
  axis tight
  yl=ylim;
  axis([xlim yl+[-1 1]*0.05*diff(yl)])

  subplot(212)
  plot(t,u2,'color',c(n,:)), hold on
  xlabel('Time [min]','fontsize',14)
  ylabel('L','fontsize',14)
  axis tight
  yl=ylim;
  axis([xlim yl+[-1 1]*0.05*diff(yl)])
