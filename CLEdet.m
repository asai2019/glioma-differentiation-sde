function xprime=CLEdet(y,params,dose)


%%%  First step estimation 
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
K6a=K6a*10;
n=10;
    
xprime=[d1+V1*CT.^3/(K1^3+CT.^3)-d1*y(1) %%%   %% PKA    % use CT stimululi substitute the role of cAMP
    V2*y(1)/(K2+y(1))-d2*y(2)   %% CREB
    1/(1+y(1)/K3)-d3*y(3)   %% PI3K
    V4*y(3)/(K4+y(3))-d4*y(4)   %% AKT
    V5*y(4)/(K5+y(4))-d5*y(5)   %% pGSK3b  % 1+SB is assumed according experiment   %(1/(1+SB/5))*(1/(1+LiCl/5))*
    (V6a)*y(6)^n/(K6a^n+y(6)^n)-d6*(1-y(5))/(K6b+(1-y(5)))*y(6) %alfa-d6*(1-y(5))^6/(K6b^6+(1-y(5))^6)   %% CyclinD1  %alfa-  V6a*y(6)/(K6a+y(6))
    V7*(y(1))/(K7+(y(1)))-d7*y(7)         %% IL6 
    V8*y(7)/(K8+y(7))-d8*y(8)   %% JAK2
    V9*y(8)/(K9^3+y(8))-d9*y(9)     %% STAT3
    V10a*y(2)/(K10a+y(2))*(y(9)/(K10c+y(9))+y(6)/(K10d+y(6)))+V10b*((1.25-y(6))/1.25*(1-y(5)))^8/(K10b^8+(1-y(5))^8)-d10*y(10)   %(TotalGSK3b-y(5)-T_GSK3b>0)* %% GFAP
    V11a*y(6)/(K11a+y(6))+V11b*y(9)/(K11b+y(9))-d11*y(11)  %% PCNA
   ];

return