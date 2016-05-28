E=2e7;
l=100;
A=10;
Iz=25;
P=10000;

Ke=[E*A/l,0,0,-E*A/l,0,0;
    0,12*E*Iz/l^3,6*E*Iz/l^2,0,-12*E*Iz/l^3,6*E*Iz/l^2;
    0,6*E*Iz/l^2,4*E*Iz/l,0,-6*E*Iz/l^2,2*E*Iz/l;
    -E*A/l,0,0,E*A/l,0,0;
    0,-12*E*Iz/l^3,-6*E*Iz/l^2,0,12*E*Iz/l^3,-6*E*Iz/l^2;
    0,6*E*Iz/l^2,2*E*Iz/l,0,-6*E*Iz/l^2,4*E*Iz/l];

Ke_local=cat(3,zeros(6),zeros(6),zeros(6));
for k=1:3
    Ke_local(:,:,k)=Ke;
end

Lambada(:,:,1)=[0,1,0,0,0,0;
                -1,0,0,0,0,0;
                0,0,1,0,0,0;
                0,0,0,0,1,0;
                0,0,0,-1,0,0;
                0,0,0,0,0,1];
            
Lambada(:,:,2)=[1,0,0,0,0,0;
                0,1,0,0,0,0;
                0,0,1,0,0,0;
                0,0,0,1,0,0;
                0,0,0,0,1,0;
                0,0,0,0,0,1];
            
Lambada(:,:,3)=[0,-1,0,0,0,0;
                1,0,0,0,0,0;
                0,0,1,0,0,0;
                0,0,0,0,-1,0;
                0,0,0,1,0,0;
                0,0,0,0,0,1];

Ke_global=cat(3,zeros(6),zeros(6),zeros(6));
for k=1:3
   Ke_global(:,:,k)=Lambada(:,:,k)'*Ke_local(:,:,k)*Lambada(:,:,k);
end

K=zeros(12);
for elem=1:3
    for r=1+(elem-1)*3:6+(elem-1)*3
        for c=1+(elem-1)*3:6+(elem-1)*3
            K(r,c)=K(r,c)+Ke_global(r-(elem-1)*3,c-(elem-1)*3,elem);
        end
    end
end

load=[0,0,0,P,0,0,0,0,0,0,0,0]';
delta=zeros(12,1);

delta(4:9,1:1)=linsolve(K(4:9,4:9),load(4:9,1:1));
load=K*delta;
for t=5:9
    load(t)=0;
end







                
