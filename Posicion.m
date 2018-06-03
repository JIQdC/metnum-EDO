source Scripts.m
global j;

lsode_options ("integration method","adams");

u0=[0 1 0 -1 1 0 -1 0];
ti=0;
tf=2*pi;

tiempos=zeros(50,5);
errores=zeros(50,5);
evf=zeros(50,5);

for i=1:100
  n=50*i+150;
  tiempos(i,1)=n;
  errores(i,1)=n;
  evf(i,1)=n;
  
  t=linspace(ti,tf,n);

  #Euler Explicito  
  j=0;
  clock=time();
  u=Theta("AtrGrav",0,u0,t);
  path=strjoin({'EE',dec2base(n,10),'.csv'},'');
  csvwrite(path,u)
  tiempos(i,2)=time()-clock;
  errores(i,2)=sqrt((1-u(n,5))^2+u(n,6)^2);
  evf(i,2)=j;
  
  #Crank-Nicolson
  j=0;
  clock=time();
  u=Theta("AtrGrav",0.5,u0,t);
  path=strjoin({'CN',dec2base(n,10),'.csv'},'');
  csvwrite(path,u);
  tiempos(i,3)=time()-clock;
  errores(i,3)=sqrt((1-u(n,5))^2+u(n,6)^2);
  evf(i,3)=j;
  
  #Euler Implicito
  j=0;
  clock=time();
  u=Theta("AtrGrav",1,u0,t);
  path=strjoin({'EI',dec2base(n,10),'.csv'},'');
  csvwrite(path,u)
  tiempos(i,4)=time()-clock;
  errores(i,4)=sqrt((1-u(n,5))^2+u(n,6)^2);
  evf(i,4)=j;
  
  #lsode
  j=0;
  clock=time();
  u=lsode("AtrGrav",u0,t);
  path=strjoin({'LS',dec2base(n,10),'.csv'},'');
  csvwrite(path,u)
  tiempos(i,5)=time()-clock;
  errores(i,5)=sqrt((1-u(n,5))^2+u(n,6)^2);
  evf(i,5)=j;
  
endfor

csvwrite("Tiempos.csv",tiempos);
csvwrite("Errores.csv",errores);
csvwrite("EvFunc.csv",evf);