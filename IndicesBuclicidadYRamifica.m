Ingreso=table2array(readtable('Sujeto 2.txt'));
Buclicidad = readfis('LDBuclicidad')
load EstadosDigra;
load RelaDigra;
load NodosDigra;

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

[Fila,~]=size(Ingreso);
Ingreso(Fila+1,1)=NaN;
JugadaReali=0;
for i=1:(Fila+1)
    if Ingreso(i,1)<10^9
        JugadaReali=JugadaReali+1;
    end
end
JugadaReali=Fila+1-JugadaReali;

Traza=[];
Indices=[];
Datosgene=[];
f=1;
k=1;

for i=1:(Fila+1)
    if Ingreso(i,1)<10^9
        for n=1:9
            Traza(k,n)=floor(Ingreso(i,1)/10^(9-n));
            Traza(k,n)=Traza(k,n)-floor(Traza(k,n)/10)*10;
        end
    else
        [Bucles,PaSinBu,PaBu,cantibucle,Fll]=bucle(Traza,EstadosDigra,251,H,Buclicidad);
        [Ramificacion,Recorinodos]=Ramifica(Traza,EstadosDigra,H1,251,H);
        [TrazaF,~]=size(Traza);
        
        if TrazaF<2
            fprintf('---NO SE HA REALIZADO MOVIMIENTOS---\n');
            return;
        end
        if Traza(1,:)~=EstadosDigra(23,:)
            fprintf('---SE DEBE EMPEZAR CON LOS CABALLOS EN EL ORDEN INICIAL---\n');
            Traza(1,:)
            fprintf('PARTIDA \n');
                f
            return;
        end
        for j=1:(length(Recorinodos)-1)
            if(mod(sum(abs(Traza(j,:)-Traza(j+1,:)))/2,2))~=((1+(-1)^j)/2)
                fprintf('---ERROR AL INTARCALAR COLOR EN LAS FICHAS---\n');
                fprintf('PARTIDA \n');
                f
                fprintf('MOVIMIENTO \n');
                j
                fprintf('Fila del TXT \n');
                i-1
                return;
            end
        end
        Indices(f,1)=floor(Bucles*1000)/1000;
        Indices(f,2)=floor(Ramificacion*1000)/1000;
        Datosgene(f,1)=PaSinBu;
        Datosgene(f,2)=PaBu;
        Datosgene(f,3)=cantibucle;
        Datosgene(f,4)=Fll;
        
        highlight(H1,Recorinodos,'EdgeColor',[1 ((1/(JugadaReali))*(f-1)) 0],"LineWidth",((1/(JugadaReali))*(f-1))*2+2);
        figure;
        plot(H,'EdgeColor','g','Layout','force','UseGravity',true);
        highlight(plot(H,'EdgeColor','g','Layout','force','UseGravity',true),Recorinodos,'EdgeColor',[1 ((1/(JugadaReali))*(f-1)) 0],"LineWidth",((1/(JugadaReali))*(f-1))*2+2)
        a=string(f);
        text(-15,15,"Partida");
        text(-11,15,a);
        f=f+1;
        Traza=[];
        k=0;
    end
    k=k+1;
end

if JugadaReali==1
    Datosgene=[Datosgene;[NaN NaN NaN NaN]];
    Indices=[Indices;[NaN NaN]];
end

figure;
b1=bar(Datosgene)
legend('Jugadas realizadas','Jugadas dentro de un bucle','Cantidad de bucles','Mínima jugadas por llegar');
xtips1 = b1(1).XEndPoints;
ytips1 = b1(1).YEndPoints;
labels1 = string(b1(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips1 = b1(2).XEndPoints;
ytips1 = b1(2).YEndPoints;
labels1 = string(b1(2).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips1 = b1(3).XEndPoints;
ytips1 = b1(3).YEndPoints;
labels1 = string(b1(3).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips1 = b1(4).XEndPoints;
ytips1 = b1(4).YEndPoints;
labels1 = string(b1(4).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
title('DATOS GENERALES')
xlabel('PARTIDAS JUGADAS')
ylabel('CANTIDAD')


figure;
b=bar(Indices);
legend('Bucles','Ramifica');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
title('INDICES BUCLICIDAD Y RAMIFICACION')
xlabel('PARTIDAS JUGADAS')
ylabel('PONDERACION')


%% Funcion Ramificación
function [Rami,ReNo] = Ramifica(Reco,EstadoT,plotnodos,Target,grafo)
RecoNodos=[];
[RecoF,~]=size(Reco);
[EstadoTF,~]=size(EstadoT);
Rami=0;
for i=1:RecoF
    for j=1:EstadoTF
        if (Reco(i,:)==EstadoT(j,:))
            RecoNodos(1,i)=j;
        end
    end
end
for i=1:RecoF-1
    RecoNodos(1,i)
    [~,numera]=size(shortestpath(grafo,RecoNodos(1,i),Target));
    numera=numera-1
    [~,denomi]=size(shortestpath(grafo,RecoNodos(1,i+1),Target))
    if numera==denomi
        Califi=1
    else
        Califi=0
    end
    %floor(numera/denomi)
    %Rami=Rami+floor(numera/denomi);
    Rami=Rami+Califi;
    fprintf('------------------------------------');
end
ReNo=RecoNodos;
Rami=Rami/(RecoF-1);

end


%% Funcion Bucle
function [recorrido,PaSinBu,PaBu,cantibucle,Fll] = bucle(Matriz3,Matriz4,Target,grafo,fuzzy)
recorrido=[];
    [f,~]=size(Matriz3);
    [k,~]=size(Matriz4);
    for  i=1:f
        for j=1:k
            if Matriz3(i,:)==Matriz4(j,:)
              recorrido(1,i)=j;
            else
            end
        end
    end
    cantibucle=0;
    Totalbucles=0;
    i=1;
    while i<=f-1
        for j=i+1:f
               if recorrido(1,i)==recorrido(1,j)
                  cantibucle=cantibucle+1;
                  fprintf('SI hay bucle en el nodo');
                  nodo=recorrido(1,i)
                  fprintf('longitud bucle');
                  longitud=j-i
                  i=j;
                  Totalbucles=Totalbucles+longitud;
                  fprintf('------------------------------------');
                  break
               end
        end
        i=i+1;
    end
    
    if cantibucle==0
        fprintf('NO hay bucle en el nodo');
    else
        cantibucle
        Totalbucles
    end
    [~,fll]=size(shortestpath(grafo,recorrido(1,f),Target));
    if fll==0
       fll=Inf; 
    else
       fll=fll-1;
    end
    
    if fll<=16
        PorcenFll=fll/16
    else
        PorcenFll=1
    end
    
    
    (Totalbucles/(f-1))
    if ((Totalbucles)==0)
        recorrido=0;
    else
        recorrido=evalfis (fuzzy, [(Totalbucles/(f-1)) PorcenFll])
    end
    PaSinBu=(f-1);
    PaBu=Totalbucles;
    Fll=fll;
end

