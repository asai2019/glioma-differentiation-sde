function xprime=CLEsto(y,params,dose,ii,GWnoise)




V1=params(1);
K1=params(2);
V2=params(3);
K2=params(4);
K3=params(5);
V4=params(6);
K4=params(7);
V5=params(8);
K5=params(9);
V6a=params(10);
K6a=params(11);
K6b=params(12);

d1=params(29);
d2=params(30);
d3=params(31);
d4=params(32);
d5=params(33);
d6=params(34);   %% CyclinD1  


    
%%%  Second step estimation 
V7=params(13);
K7=params(14);
V8=params(15);
K8=params(16);
V9=params(17);
K9=params(18);

V10a=params(19);
K10a=params(20);
V10b=params(21);
K10b=params(22);

K10c=params(23);
K10d=params(24);

V11a=params(25);
K11a=params(26);
K11b=params(27);
V11b=params(28);


d7=params(35);
d8=params(36);
d9=params(37);
d10=params(38);   %% GFAP
d11=params(39);   %% PCNA

CT = dose;
% K6a=K6a*10;
TotalGSK3b=1;  % Threshold of GSK3beta
    
xprime=[sqrt(abs(d1+V1*CT.^3/(K1^3+CT.^3)))*GWnoise(ii,1)-sqrt(abs(d1*y(1))*GWnoise(ii,2)) %%%   %% PKA    % use CT stimululi substitute the role of cAMP
    sqrt(abs(V2*y(1)/(K2+y(1)))*GWnoise(ii,3))-sqrt(abs(d2*y(2))*GWnoise(ii,4))   %% CREB
    sqrt(abs(1/(1+y(1)/K3))*GWnoise(ii,5))-sqrt(abs(d3*y(3))*GWnoise(ii,6))   %% PI3K
    sqrt(abs(V4*y(3)/(K4+y(3)))*GWnoise(ii,7))-sqrt(abs(d4*y(4))*GWnoise(ii,8))   %% AKT
    sqrt(abs(V5*y(4)/(K5+y(4)))*GWnoise(ii,9))-sqrt(abs(d5*y(5))*GWnoise(ii,10))   %% pGSK3b  % 1+SB is assumed according experiment   %(1/(1+SB/5))*(1/(1+LiCl/5))*
    sqrt(abs((V6a+0.005)*y(6)^10/(K6a^10+y(6)^10))*GWnoise(ii,11))-sqrt(abs(d6*(TotalGSK3b-y(5))/(K6b+(TotalGSK3b-y(5)))*y(6))*GWnoise(ii,12)) %alfa-d6*(1-y(5))^6/(K6b^6+(1-y(5))^6)   %% CyclinD1  %alfa-  V6a*y(6)/(K6a+y(6))
    sqrt(abs(V7*(y(1))/(K7+(y(1))))*GWnoise(ii,13))-sqrt(abs(abs(d7*y(7))*GWnoise(ii,14)))         %% IL6 
    sqrt(abs(V8*y(7)/(K8+y(7)))*GWnoise(ii,15))-sqrt(abs(d8*y(8))*GWnoise(ii,16))   %% JAK2
    sqrt(abs(V9*y(8)/(K9^3+y(8)))*GWnoise(ii,17))-sqrt(abs(d9*y(9))*GWnoise(ii,18))     %% STAT3
    sqrt(abs(V10a*y(2)/(K10a+y(2))*(y(9)/(K10c+y(9))+y(6)/(K10d+y(6)))+V10b*((1.25-y(6))/1.25*(TotalGSK3b-y(5)))^8/(K10b^8+(TotalGSK3b-y(5))^8))*GWnoise(ii,19))-sqrt(abs(d10*y(10))*GWnoise(ii,20))   %(TotalGSK3b-y(5)-T_GSK3b>0)* %% GFAP
    sqrt(abs(V11a*y(6)/(K11a+y(6))+V11b*y(9)/(K11b+y(9)))*GWnoise(ii,21))-sqrt(abs(d11*y(11))*GWnoise(ii,22))  %% PCNA
   ];

return