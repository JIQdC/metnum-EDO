source Scripts.m
global j;

lsode_options ("integration method","adams");

u0=[0 0.5 0 -0.5 1 0 -1 0];
ti=0;
tf=2.71414;

for i=[10,20,50,75,100,200,300]
  n=50*i+150;
  
  t=linspace(ti,tf,n);
  
  e=zeros(n,5);
  e(:,1)=t;

  #Euler Explicito  
  u=Theta("AtrGrav",0,u0,t);
  e(:,2)=EMec(u,t);

  #Crank-Nicolson
  u=Theta("AtrGrav",0.5,u0,t);
  e(:,3)=EMec(u,t);
  
  #Euler Implicito
  u=Theta("AtrGrav",1,u0,t);
  e(:,4)=EMec(u,t);
  
  #lsode
  u=lsode("AtrGrav",u0,t);
  e(:,5)=EMec(u,t);
  
  path=strjoin({'En',dec2base(n,10),'.csv'},'');
  csvwrite(path,e)
endfor