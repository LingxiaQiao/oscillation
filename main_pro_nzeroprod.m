function [period,xmax,xmin, ss]=main_pro_nzeroprod(matrix_v,matrix_K,vector_r,vector_delta,Hill_n,Jplus,Jabs,starting_point)
% calcualte deterministic behavior

% Input variables:
%     matrix_v,matrix_K,vector_r,vector_delta are kinetic parameters
%     Hill_n is Hill coefficient 
%     Jplus,Jabs are network topologies information
%     starting_point is set to be 2, i.e, we selet a whole period when node B reaches maximal value

% Output variables:
%     period:period length 
%     xmax: maximal values of nodes A, B, C in a complete period
%     xmin: minimal values of nodes A, B, C in a complete period
%     ss: the state in a complete period that node B reaches maximal value


period=0;xmax=[0 0 0]; xmin=[0 0 0];ss=[0 0 0 ];
Nnode=size(matrix_v,1);
T=1000;

[TT,Y]=ode15s(@(t,x)myfun(t,x,matrix_v,matrix_K,vector_r,vector_delta,Hill_n,Jplus,Jabs),[0 T],zeros(1,Nnode));% figure;plot(TT,Y)
in=find(abs(TT-0.7*T)<0.1*T); 
if isempty(in)  
    return;
end
s=0;
for i=1:length(matrix_v)
    if max(abs(Y(end,i)-Y(in(1):end,i))-Y(end,i)/10)<1e-5
       s=s+1; 
    end
end
if s==length(matrix_v)
    return;
end
T_detail=TT(floor(length(Y)/2):end);
Y_detail=Y(floor(length(Y)/2):end,starting_point);Y_detail_all=Y(floor(length(Y)/2):end,:);
[pks,locs]=findpeaks(Y_detail);


if length(locs)<=2
    return;
end
per=diff(T_detail(locs));
if (max(per)-min(per)<min(per)/10)&&(max(pks)-min(pks)<min(pks)/10)
    period=mean(per);
    amplitude=pks(1)-min(Y_detail);
    ts=T_detail(locs(1):locs(2));
    xs=Y_detail_all(locs(1):locs(2),:); 
    xmax=max(xs);  xmin=min(xs); ss=xs(1,:);
end


end


function dx=myfun(t,x,matrix_v,matrix_K,vector_r,vector_delta,Hill_n,Jplus,Jabs)
Nnode=size(matrix_v,1);
cc=ones(Nnode,1)*x';
a=1+sum(Jabs.*(cc./matrix_K).^Hill_n,2);%vector
b=  sum(Jplus.*matrix_v.*(cc./matrix_K).^Hill_n,2); %vector
dx=-vector_r.*x+(b+vector_delta)./a+0.01;
end


