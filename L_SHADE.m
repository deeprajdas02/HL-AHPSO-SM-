function [xo,fo,FEs,X,Y,phase] = L_SHADE(dim,X,Fx,Tmax,MaxFEs,FEs,N,xmin,xmax,H,Co,Fo,phase,xo,fo, fhd, func_no, C)
% Initialize dynamic parameters
disp("LShade");
Ni=N; 
MCR=Co*ones(1,H);
MF=Fo*ones(1,H);
km=1;
Fi=zeros(N,1);
Cr=Fi;
dF=Fi;
t=0;
NA=round(2.6*N); 
Y=zeros(NA,dim);
ny=0;

% Initialize archive
FEs = FEs + N; 
k=0;
fo_prev = inf;  
flag = 0; 
u=zeros(1,dim);

while FEs < MaxFEs
    k=k+1;
    nb=max(2,ceil(0.11*N));
    for i = 1:N
        j=randi(H);
        Fi(i)=0;
        while Fi(i)<=0
            Fi(i)=min(1,MF(j)+0.1*tan(pi*(rand-0.5)));
        end
        rb=randperm(nb,2);
        if rb(1)==i, rb(1)=rb(2); end 
        r1=randperm(N,2);
        if r1(1)==i, r1(1)=r1(2); end 
        r2=randperm(N+ny,2);
        if r2(1)==i, r2(1)=r2(2); end 
        if r2(1)<=N
            V(i,:)=X(i,:)+Fi(i)*(X(rb(1),:)-X(i,:)) ...
                +Fi(i)*(X(r1(1),:)-X(r2(1),:));
        else
            V(i,:)=X(i,:)+Fi(i)*(X(rb(1),:)-X(i,:)) ...
                +Fi(i)*(X(r1(1),:)-Y(r2(1)-N,:));
        end  
        V(i,:) = max(min(V(i,:),xmax),xmin);
    end
  
    for i = 1:N
        j=randi(H);
        Cr(i)=max(0,min(1,MCR(j)+0.1*randn));
        di = randi(dim);
        for j = 1:dim
            if (j == di) || (rand <= Cr(i))         
                u(1,j)=V(i,j);
            else
                u(1,j)=X(i,j);
            end
            u = max(min(u,xmax),xmin);
        end
        Fu = feval(fhd,u',func_no,C);
        if Fu < Fx(i)
            if ny<NA 
                ny=ny+1;
                j=ny;
            else
                j=randi(NA);
            end
            Y(j,:)=X(i,:);
            dF(i)=Fx(i)-Fu;
            X(i,:)=u;
            Fx(i)=Fu;
        else
            Cr(i)=-1;
        end
     end
    [Fx, ids]=sort(Fx);
    X=X(ids,:);
    if fo > Fx(1)
        fo=Fx(1); 
        xo = X(1,:);
    end

    FEs = FEs + N;
    t=t+1;
    Co1=0;
    Co2=0;
    Fo1=0;
    Fo2=0;
    sdF=0; 
    for i=1:N
      if Cr(i)>= 0
        sdF=sdF+dF(i);
      end
    end
    j=0;
    for i=1:N
        if Cr(i)>= 0
            j=j+1;
            di=dF(i)/sdF;
            Co1=Co1+di*Cr(i);
            Co2=Co+di*Cr(i)^2;
            Fo1=Fo1+di*Fi(i);
            Fo2=Fo2+di*Fi(i)^2;
         end
    end
    if j>0
        Co=Co2/Co1; 
        Fo=Fo2/Fo1;
    end
    MF(km)=Fo;
    MCR(km)=Co;
    km=km+1;
    if km>H
        km=1;
    end
    N=round((4-Ni)/MaxFEs*FEs+Ni);
    NA=round((1+(MaxFEs-FEs)/MaxFEs*1.6)*N);
    if ny > NA
       ny=NA;
    end
    
    if fo_prev == fo 
        flag = flag + 1;
    else 
        flag = 0;
        fo_prev = fo;
    end
    if flag == 20 || t == Tmax
        phase = "AHPSO";
        disp("returned");
        return
    end
end

return

