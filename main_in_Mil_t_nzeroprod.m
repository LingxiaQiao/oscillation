function [tau1,tau2 ,tau3]=main_in_Mil_t_nzeroprod(matrix_v,matrix_K,vector_r,vector_delta,Hill_n,Jplus,Jabs,upper_time,ss,n,period)
% intrinisc noise is considered and calcualte dimensionless autocorrelation time, i.e., ratios of autocorrelation time to period 

% Input variables:
%     matrix_v,matrix_K,vector_r,vector_delta are kinetic parameters
%     Hill_n is Hill coefficient 
%     Jplus,Jabs are network topologies information
%     upper_time is the ending point of the time variable
%     ss is the initial state of nodes A, B, and C
%     n is the number of iterations
%     period is the period length

% Ouptu variables: tau1, tau2 and tau3 are dimensionless autocorrelation time for nodes A,B
% and C, respectively.

tau1=0;tau2=0;tau3=0;
Nnode=size(matrix_v,1);
T=ceil(upper_time);dt=T/n;
rand_num=normrnd(0,1,n,Nnode)*sqrt(dt);
Y=ones(n,Nnode);
V=100;   Y(1,:)=V*ss;  % Y is the noisy trajectory 
% simulate the dynamics
for j=2:n+1
    x=Y(j-1,:)';
    cc=ones(Nnode,1)*x';
    a=V^3+sum(Jabs.*(cc./matrix_K).^Hill_n,2);
    b=sum(Jplus.*matrix_v.*(cc./matrix_K).^Hill_n,2);  % vector
      
    c=[0.01*ones(Nnode,1)*V,  V*(b+V^3*vector_delta)./a,    -vector_r.*x]; 
    sigma_Xn=sqrt(sum(abs(c),2));
    
    temp=Y(j-1,:)'+sum(c,2)*dt+sigma_Xn.*rand_num(j-1,:)'+...
        0.5*sigma_Xn.^2.*((rand_num(j-1,:)').^2-dt);
    
    Y(j,:)=(abs(temp)+temp)/2;  
end
    
% calculate dimensionless autocorrelation time
Y=Y/V;
t=(0:dt:T)';
for indd=1:3
    
     d_ind=ceil(period/(t(2)-t(1)));
     temp=Y((end-d_ind):end,indd);
     me=mean(Y(:,indd));
     acf=mean((Y(1:(end-d_ind),indd)-me).*(Y((1+d_ind):end,indd)-me));
     s=mean((Y(1:end,indd)-me).*(Y(1:end,indd)-me));
     
     d_ind=ceil(period/(t(2)-t(1)));
     if ~isempty(find(autocorr(Y(:,indd),'NumLags',d_ind-1)<0))
         if acf>0
             tau=-1/log((acf/s));
         else
             tau=0;
         end
     else
         tau=0;
     end
     
     
    eval(['tau',num2str(indd),'=tau;']);
    eval(['acf',num2str(indd),'=acf;']);
end


% plot autocorrelation function and noisy trajectory
for indd=1:3
    d_ind=ceil(3*period/(t(2)-t(1)));
    acf=autocorr(Y(:,indd),'NumLags',d_ind-1);
    acf=acf.*length(Y)./(length(Y)-(0:(d_ind-1)))';
    t_cal=t(1:d_ind);
    warning('off');[acf_peak,loc]=findpeaks(acf,t_cal,'MinPeakHeight',0,'MinPeakDistance',period/2);
    eval(['acf',num2str(indd),'=acf;']);
end


figure;set(gcf,'unit','centimeters','position',[2,2,30,12]);

subplot(1,2,1);hold on;  % autocorrelation function
plot(t_cal(10:30:end),exp(-t_cal(10:30:end)./(tau2*period)).*cos(2*pi*t_cal(10:30:end)/period),'og','linewidth',0.5);
h2=plot(t_cal,acf2,'g','linewidth',0.5);
legend([h2],'B');
ylim([-1 1])

subplot(1,2,2);plot(t,Y(:,2));title('B');  % noisy trajectory


end

