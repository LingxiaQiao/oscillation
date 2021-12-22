function netplot(A)
% plot a network topology
 set(0,'DefaultLineLineWidth',3);
% 
% %(axes &legend) font    ,label
 set(0,'DefaultAxesFontSize',20,'DefaultAxesFontWeight','bold');
% 
% %text
 set(0,'DefaultTextFontSize',20,'DefaultTextFontWeight','bold');

if (size(A,1)==1&size(A,2)==9)
    A=reshape(A,3,3);
end

      
x=[0  1   1  0];
y=[1  1   0  0];  
%plot(x,y,'r*');
%

axis([-0.2 1 -0.2 1]);
title('topology'); 
text([0.05],[0.95],'A')
%text([0.85],[1.2],'B')
text([0.05],[0.05],'B')
text([0.72],[0.05],'C');
%%  A->       B->        C->      D->
h_axes = get(gcf,'CurrentAxes');    %get axes handle.
axesoffsets = get(h_axes,'Position');

% A->A (start,end)    ;B->A (start,end)   ; C->A (start,end)  X1
% A->B (start,end)    ;B->B(start,end)    ; C->B (start,end)  X2
% A->C (start,end)    ;B->C (start,end)   ; C->C (start,end)  X3

rr=0.6;rl=0.3;
X1=[0.2  0.28;0.30 0.30;0.62 0.4 ;]/rr-rl;
Y1=[0.8   0.7;0.35 0.65;0.40 0.65;]/rr-rl;

X3=[0.37 0.60;0.40 0.60;0.8 0.7;]/rr-rl;
Y3=[0.62 0.38;0.35 0.35;0.15 0.25;]/rr-rl;

X2=[0.35 0.35;0.2  0.28;0.60 0.40;]/rr-rl;
Y2=[0.65 0.35;0.15 0.25;0.30 0.30;]/rr-rl; 


X1=X1*axesoffsets(3)+axesoffsets(1); 
X2=X2*axesoffsets(3)+axesoffsets(1); 
X3=X3*axesoffsets(3)+axesoffsets(1); 

Y1=Y1*axesoffsets(4)+axesoffsets(2);
Y2=Y2*axesoffsets(4)+axesoffsets(2); 
Y3=Y3*axesoffsets(4)+axesoffsets(2);

X=zeros(3,2,3);
X(:,:,1)=X1;X(:,:,2)=X2;X(:,:,3)=X3;

Y=zeros(3,2,3);
Y(:,:,1)=Y1;Y(:,:,2)=Y2;Y(:,:,3)=Y3;
for i=1:3
    for j=1:3
        if A(i,j)>0
           annotation('arrow', X(j,:,i), Y(j,:,i)); 
        end
        if A(i,j)<0
           annotation('arrow',X(j,:,i), Y(j,:,i),'HeadStyle','ellipse');
         end
    end
end
axis off;
set(0,'DefaultLineLineWidth','remove');
% 
% %(axes &legend) font    ,label
 set(0,'DefaultAxesFontSize','remove','DefaultAxesFontWeight','remove');
% 
% %text
 set(0,'DefaultTextFontSize','remove','DefaultTextFontWeight','remove');
