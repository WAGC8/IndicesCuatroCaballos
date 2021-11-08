%---------------------Genera matriz Estados con todas las posibidades de
%ordenar las fichas en el tablero 8x7x6x5=1680--------------------------
k=1;
Estados=zeros(1680,9);
for i=1:8
    for j=1:210
            Estados(k,i)=1;
            k=k+1;
    end
end
k=1;
while k<=1680
    for i=1:8
        if (k~=1681)&&(Estados(k,i)==0)
            for j=1:30
                Estados(k,i)=2;
                k=k+1;
            end
        end
    end
end
k=1;
while k<=1680
    for i=1:8
        if (k~=1681)&&(Estados(k,i)==0)
            for j=1:5
                Estados(k,i)=3;
                k=k+1;
            end
        end
    end
end
k=1;
while k<=1680
    for i=1:8
        if (k~=1681)&&(Estados(k,i)==0)
                Estados(k,i)=4;
                k=k+1;
        end
    end
end
for i=1:1680
    Estados(i,9)=Estados(i,8);
    Estados(i,8)=Estados(i,7);
    Estados(i,7)=Estados(i,6);
    Estados(i,6)=Estados(i,5);
    Estados(i,5)=0;
end

%---Genera matriz Rela(Matriz de adyacencia de los 1680 nodos)-------------
%Rela matriz simetrica genera grafo no diriguido
Rela=zeros(1680);
Color=zeros(1680);
Posi=zeros(2);
for i=1:1680
    for j=1:1680
        Matrizi=ConvertirMatriz(Estados,i);
        Matrizj=ConvertirMatriz(Estados,j);
        cont=0;
        for m=1:3
            for n=1:3
                if Matrizi(m,n)~=Matrizj(m,n)
                    cont=cont+1;
                    FiJu=abs(Matrizi(m,n)+Matrizj(m,n));
                end
            end
        end
        
        if cont==2
            h=1;
            multicero=0;
            RestaMatriz=Matrizi-Matrizj;
            for m=1:3
                for n=1:3
                    if Matrizi(m,n)*Matrizj(m,n)~=0
                        multicero=multicero+1;
                    end
                    
                    if RestaMatriz(m,n)~=0
                        Posi(h,1)=m;
                        Posi(h,2)=n;
                        h=h+1;
                    end
                end
            end
            
            if ((Posi(1,1)-Posi(2,1))^2+(Posi(1,2)-Posi(2,2))^2==5)&&(multicero==3)
               Rela(i,j)=1; 
               if mod(FiJu,2)==0
                   Color(i,j)=FiJu;
               else
                   Color(i,j)=FiJu;
               end
            end
            
        end
    end
end

G=graph(Rela);

tamano=1680;
TurnoBN=zeros(tamano,1);
for i=1:tamano
  TurnoBN(i,1)=-1;
end
TurnoBN(129,1)=0;

for k=1:19
    for i=1:tamano
        if (TurnoBN(i,1)==(k-1))
            for j=1:tamano
                if (TurnoBN(j,1)==-1)&&(Color(i,j)~=0)&&(mod(Color(i,j),2)==((1+(-1)^k)/2))
                    TurnoBN(j,1)=k;
                end
            end
        end
    end
end

contador=0;
for i=1:tamano
     if TurnoBN(i,1)~=-1
     contador=contador+1;
     end
end

RelaDigra=Rela;

for i=1:1680
    if (TurnoBN(i,1)~=-1)
        for j=1:1680
            if (RelaDigra(i,j)~=0)&&(mod(Color(i,j),2)~=mod(TurnoBN(i,1),2))
              RelaDigra(i,j)=0;
            end
        end
    end
end
RelaDigra(115,:)=0;
RelaDigra(418,:)=0;
RelaDigra(591,:)=0;
RelaDigra(702,:)=0;
RelaDigra(985,:)=0;
RelaDigra(1091,:)=0;
RelaDigra(1269,:)=0;
RelaDigra(1572,:)=0;

NodosDigra=[];
EstadosDigra=[];
posi=0;
for i=1:1680
    if (TurnoBN(i,1)~=-1)
        posi=posi+1;
        NodosDigra(posi)=i;
        EstadosDigra(posi,:)=Estados(i,:);
    end
end

G1=digraph(RelaDigra);
H = subgraph(G1,NodosDigra);
H1 = plot(H,'EdgeColor','g','Layout','force','UseGravity',true);
labelnode(H1,23,'Start')
labelnode(H1,251,'End')
highlight(H1,[23 251],'NodeColor','g')
pathway = shortestpath(H,23,251);
pathway1 = shortestpath(H,251,23);
highlight(H1,pathway,'EdgeColor','cyan',"LineWidth",1.5)
highlight(H1,pathway1,'EdgeColor','cyan',"LineWidth",1.5)
set(gca,'color',[0.066 0.066 0.066])
colormap jet

%% Funcion CONVERTIR vector a matriz

function Tablero = ConvertirMatriz(vector,n)
    Tablero=[];
    for i=1:3
        for j=1:3
          Tablero(i,j)=vector(n,j);
          if i==2
              Tablero(i,j)=vector(n,j+3);
          end  
          if i==3
              Tablero(i,j)=vector(n,j+6);
          end
        end
    end
end




