1;

function xf=AtrGrav(x,t)
  global j;
  j++;
  m1=4;
  m2=4;
  nm=sqrt((x(7)-x(5))^2+(x(8)-x(6))^2)**3;
  xf=zeros(8,1);
  xf(5:8)=x(1:4);
  xf(1)=m2*(x(7)-x(5))/nm;
  xf(2)=m2*(x(8)-x(6))/nm;
  xf(3)=m1*(x(5)-x(7))/nm;
  xf(4)=m1*(x(6)-x(8))/nm;
endfunction

function u=Theta(f,theta,u0,t)
  #Aplica el metodo theta para resolver un sistema de ec diferenciales.
  #f es la funcion caracteristica del sistema
  #theta es el parametro del metodo
  #x0 es el conjunto de las condiciones iniciales (fila)
  #t es el vector de tiempo en el que se quiere obtener la solucion (columna)
  
  n=size(t)(2);   #puntos donde evaluo la solucion
  m=size(u0)(2);  #tamaño del sistema
  
  u=zeros(n,m);   #matriz solucion: en cada fila tengo el valor de las incognitas,
                  #y las columnas me marcan el paso del tiempo
  
  u(1,:)=u0;       #el primer punto es la condicion inicial
  
  for i = 2:n
    #defino paso de tiempo
    h=t(i)-t(i-1);
    
    #Euler explicito
    u(i,:)=u(i-1,:)+h*feval(f,u(i-1,:)',t(i-1))';
    
    #si theta es cero, entonces se termina la iteracion
    if theta==0
      continue
    endif
        
    #si theta no es cero, tengo un metodo implicito. Itero usando como semilla
    #el valor de Euler explicito.
    e=1;
    k=0;
    while (e>1E-12 && k<50)
      k++;
      ub=u(i,:);
      
      #calculo nuevamente u(i,:) pero ahora uso el metodo theta
      k1=feval(f,u(i-1,:)',t(i-1))';
      k2=feval(f,u(i,:)',t(i))';
      u(i,:)=u(i-1,:)+h*((1-theta)*k1+theta*k2);
      
      #saco error
      e=norm(u(i,:)-ub);
    end
  endfor
endfunction

function e=EMec(u,t)
  #calcula la energia mecanica del sistema de masas m1 y m2 del sistema.
  #se computa la energia por unidad de masa (dado que las masas son iguales)
  #dado el resultado de la resolucion del sistema de EDOs en la matriz u, computa
  #la energia para cada t
  
  n=size(t)(2);
  
  e=zeros(n,1);
  
  for i=1:n
    k=2*(u(i,1)^2+u(i,2)^2+u(i,3)^2+u(i,4)^2);   #cinetica
    p=16/norm(u(i,5:6)-u(i,7:8));                 #potencial
    e(i,1)=k-p;    #mecanica=cinetica+potencial
  endfor
endfunction   