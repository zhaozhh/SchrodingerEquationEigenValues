%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Zhao Zhen-Hua 
% Date  : 2014-6-20
% Purpose: Searching the eigenvalues of Schordinger equation.
% The Schordinger equation has the form f''+(E-U))f=0,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U=load('Vzk=1C2=20C3=0p40p50p6.txt');
 x=U(:,1);%coordinate
% x=complementX(x,1);
 y=U(:,3);%potential.
% y=complementX(y,2);
U=[x,y];
figure
plot(x,y)
title('U(z)')

eigenRange=linspace(min(y)+0.1,70,1000);
[eigenValues,waveFunctions]=eigenNumerov21(eigenRange,U);
%The first column of waveFunctions is coordinate x.

fprintf('The eigenvlue is:  %3.16f, \n', eigenValues);
fprintf('the max value of U is :  %3.16f, \n', max(y));

for k=1:length(eigenValues(:,1))
    %hold all
    figure
    plot(x,waveFunctions(:,k+1))
    title(['E=',num2str(eigenValues(k,1))])
end
 


save([pwd,'/eignvaluek=1C2=20C3=0p5.txt'],'eigenValues','-ascii','-double')
save([pwd,'/wavesk=1C2=20C3=0p5.txt'],'waveFunctions','-ascii','-double')







