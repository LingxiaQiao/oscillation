load data_in.mat

% check whether J1-J5 are the five representive networks 
figure;
subplot(1,5,1);netplot3(J_1)
subplot(1,5,2);netplot3(J_2)
subplot(1,5,3);netplot3(J_3)
subplot(1,5,4);netplot3(J_4)
subplot(1,5,5);netplot3(J_5)

set(0,'DefaultLineLineWidth',1);set(0,'DefaultAxesFontSize',28,'DefaultAxesFontWeight','bold','DefaultAxesFontName','Arial');set(0,'DefaultTextFontSize',28,'DefaultTextFontWeight','bold','DefaultTextFontName','Arial');

for cc=1:5     
    
    matrix_v=matrix_v_all(:,:,cc);     matrix_K=matrix_K_all(:,:,cc);
    vector_r=vector_r_all(:,cc);   vector_delta=vector_delta_all(:,cc);
    
    eval(['J=J_',num2str(cc),';']);
    Jabs=abs(J);Jplus=(J+Jabs)/2;Jminus=(J-Jabs)/2;
    Hill_n=3;Nnode=size(J,1);  node_nega=zeros(Nnode,1);
    
    % deterministic behavior
    [period,xa,xi,ss]=main_pro_nzeroprod(matrix_v,matrix_K,vector_r,vector_delta,Hill_n,Jplus,Jabs,2);

    % add intrinsic noise
    [tau1,tau2 ,tau3]=main_in_Mil_t_nzeroprod(matrix_v,matrix_K,vector_r,vector_delta,Hill_n,Jplus,Jabs,period*100,ss,10^5,period);
    disp(['Dimensionless autocorrelation time for nodes A, B and C:',num2str([tau1,tau2 ,tau3])]);
   
    xlim([0 25*period]);ylim([0 25])

end


            