
%------------------------------------ Bump --------------------------------------------------

%% LF
clear all
close all
clc

N_ele = [1498, 5998, 23998, 95998]';

lc = [0.02, 0.01, 0.005, 0.0025]';

N2_s = [3.28925764E-03, 2.37568771E-03, 1.54114026E-03, 9.92628629E-04]';

Results = table(N_ele,lc,N2_s,...
    'VariableNames', {'Numero elementi','Lunghezza caratteristica','Norma-2 entropia'});
display(Results)

for i=1:3
    p(i,1) = log(N2_s(i)/N2_s(i+1))/log(lc(i)/lc(i+1));
end

Order = table([12,23,34]',p,...
    'VariableNames', {'Confronto mesh','p'});
display(Order)

% Richardson effettivo

uh = N2_s(4);
u2h = N2_s(3);
u4h = N2_s(2);
p = log((u4h-u2h)/(u2h-uh))/log(2);
u0 = uh + (u2h-uh)/(2^p-1);
Eh = abs(uh - u0);
Fs = 3;
GCI = Eh*Fs;

richardson_effettivo = table(p,u0,Eh,GCI);

display(richardson_effettivo)

% Richardson teorico

p = 1;
u0 = (2*p*uh-u2h)/(2^p-1);
Eh = uh-u0;
GCI = Eh*Fs;

richardson_teorico = table(p,u0,Eh,GCI);
display(richardson_teorico)

figure(1)
plot(lc,N2_s,'-o',LineWidth=1.25);
grid on
xlabel('lc')
ylabel('Norma-2 entropia')
hold on 

%% ROE
clear all
%close all
clc

N_ele = [1498, 5998, 23998, 95998]';

lc = [0.02, 0.01, 0.005, 0.0025]';

N2_s = [2.81572551E-03,1.78647193E-03 , 1.01367244E-03, 0.50904E-03]';

Results = table(N_ele,lc,N2_s,...
    'VariableNames', {'Numero elementi','Lunghezza caratteristica','Norma-2 entropia'});
display(Results)

for i=1:2
    p(i,1) = log(N2_s(i)/N2_s(i+1))/log(lc(i)/lc(i+1));
end

Order = table([12,23]',p,...
    'VariableNames', {'Confronto mesh','p'});
display(Order)

% Richardson effettivo

uh = N2_s(3);
u2h = N2_s(2);
u4h = N2_s(1);
p = log((u4h-u2h)/(u2h-uh))/log(2);
u0 = uh + (u2h-uh)/(2^p-1);
Eh = abs(uh - u0);
Fs = 3;
GCI = Eh*Fs;

richardson_effettivo = table(p,u0,Eh,GCI);

display(richardson_effettivo)

% Richardson teorico

p = 1;
u0 = (2*p*uh-u2h)/(2^p-1);
Eh = abs(uh-u0);
GCI = Eh*Fs;

richardson_teorico = table(p,u0,Eh,GCI);
display(richardson_teorico)

figure(1)
plot(lc,N2_s,'-o',LineWidth=1.25);
xlabel('lc')
ylabel('Norma-2 entropia')
legend('LF','ROE','Location','NorthWest')

%------------------------------------ Intake ------------------------------------------------
%%





%------------------------------------ LS59 --------------------------------------------------
%%