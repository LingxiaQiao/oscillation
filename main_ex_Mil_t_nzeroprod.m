function [tau1,tau2 ,tau3]=main_ex_Mil_t_nzeroprod(matrix_v,matrix_K,vector_r,vector_delta,Hill_n,Jplus,Jabs,upper_time,ss,n,period)
% extrinisc noise is considered and calcualte dimensionless autocorrelation time, i.e., ratios of autocorrelation time to period 

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
T=upper_time;dt=T/n;
ebs=0.5;  %ebs=0.1;
rand_v=normrnd(0,1,Nnode,Nnode,n)*sqrt(dt);
rand_r=normrnd(0,1,Nnode,1,n)*sqrt(dt);
rand_d=normrnd(0,1,Nnode,1,n)*sqrt(dt);

Y=ones(n,Nnode);Y(1,:)=ss;  % Y is the noisy trajectory 

eta  =zeros(Nnode,Nnode);
eta_r=zeros(Nnode,1);
eta_d=zeros(Nnode,1);
% simulate the dynamics
for j=2:n+1
    x=Y(j-1,:)';                  
    cc=ones(Nnode,1)*x';
    a= 1 +sum(Jabs.*(cc./matrix_K).^Hill_n,2);%vector
    

    eta=   eta -dt*eta  +ebs*reshape(rand_v(:,:,j-1),3,3)+ 0.5*ebs^2*(reshape(rand_v(:,:,j-1),3,3).^2-dt);
    eta_r=eta_r-dt*eta_r+ebs*reshape(rand_r(:,:,j-1),3,1)+ 0.5*ebs^2*(reshape(rand_r(:,:,j-1),3,1).^2-dt);
    eta_d=eta_d-dt*eta_d+ebs*reshape(rand_d(:,:,j-1),3,1)+ 0.5*ebs^2*(reshape(rand_d(:,:,j-1),3,1).^2-dt);
    matrix_dx=[(Jplus.*matrix_v.*(ones(Nnode,Nnode)+eta).*(cc./matrix_K).^Hill_n)./(a*ones(1,Nnode)),vector_delta.*(1+eta_d)./a,-vector_r.*(1+eta_r).*x];
    
    temp=Y(j-1,:)'+(sum(matrix_dx,2)+0.01)*dt;
    
    Y(j,:)=(abs(temp)+temp)/2;  
end

t=(0:dt:T)';

% calculate dimensionless autocorrelation time
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

subplot(1,2,1);hold on; % autocorrelation function
plot(t_cal(10:30:end),exp(-t_cal(10:30:end)./(tau2*period)).*cos(2*pi*t_cal(10:30:end)/period),'og','linewidth',0.5);
h2=plot(t_cal,acf2,'g','linewidth',0.5);
legend([h2],'B');

subplot(1,2,2);plot(t,Y(:,2));title('B'); % noisy trajectory


end
