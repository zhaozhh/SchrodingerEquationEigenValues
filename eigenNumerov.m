function [eigenValues,waveFunctions]=eigenNumerov(eigenRange,potential)
% 2017.08.16, small changes be made from eigenNumerov2.m
% eigenRange: is a vector which includes a series of numbers in the
%             range of eigenvalues.
% potential: which is the potential of Schordinger equation and is a vector.
%           The first column of potential is coordiante x, the second
%           column is the values of potential. x runs form 0 to its end.
% eigenValues: its first column is the eigen values, its second column is
%              the even or odd property of its corresponding wave function. 
% 
% The Schordinger equation has the form f''+(E-U))f=0,
% 
% Author: Zhao Zhen-Hua,2014-6-20.
% Algorithm from Tao Pang P_107.

global U nx h

U = potential(:,2);
x = potential(:,1);
nx=numel(U);
%fprintf('nx=%d\n',nx);
h=x(2)-x(1);
eigenValues=[];
waveFunctions=[];
%options = optimset('TolX',1e-10);

f0=zeros(length(eigenRange),1);
% positionEn=[];

i=1;
for en=eigenRange
    f0(i)=fen(en);
    i=i+1;
end

for j=1:length(eigenRange)-1
    if -sign(f0(j))==sign(f0(j+1))%changed at 2017.08.16 
        en=eigenRange(j);
        [m2,xx]=fzero(@fen,en);%,options);%
        if abs(xx)<0.0001%
            [~, ~,u,~,~]=wave(m2);    %changed at 2017.08.16        
            eigenValues=[eigenValues;m2];% Saving.
            waveFunctions=[waveFunctions, u'];
            % fprintf('The eigenvlue  are:  %f \n', eigenValues)
        end
    end
    
end
%  fprintf('yr(i):  %f \n', yr)

% 
[eigenValues,ix]=sort(eigenValues(:,1));%sort
waveFunctionsSort=zeros(size(waveFunctions));

for k=1:length(ix)
    waveFunctionsSort(:,k)=waveFunctions(:,ix(k));
end

waveFunctions=[x,waveFunctionsSort];


function f0=fen(en)
global  h  nx
[ul,ur,~,im,~]=wave(en);%changed at 2017.08.16 
f0 = (ul(im+1)-ul(im-1))+(ur(nx-im+2)-ur(nx-im));
f0 = f0/(2*h*ur(nx-im+1));
%ur(nx-im+1)-ul(im)

% Method to calculate the wavefunction.
function [ul,ur,u,im,ratio]=wave(en)
global  U  h nx 
y = zeros(1,nx);
s = zeros(1,nx);
u = zeros(1,nx);
u0 = 0; 
u1 = 0.01;
%Set up function q(x) in the equation
ql = (en-U);
qr = fliplr(ql);
%fprintf('length(ql)= %d \n', length(ql));
%fprintf('length(U)= %d \n', length(U));

% Find the matching point at the right turning point
im=2;%if the value of im not changed,there are not matching point
for i=1:nx-1  
    
    if (-sign(ql(i))==sign(ql(i+1))) && ql(i)>0%changed at 2017.08.16 
        im = i;
        break;
    end
    
end
%fprintf('im=%d\n',im)
% Carry out the Numerov integrations

nl = im+1 ;
nr = nx-im + 2;
ul = numerov(nl, h, u0, u1, ql, s);
ur = numerov(nr, h, u0, u1, qr, s);

% Find the wavefunction on the left
ratio = ur(nx-im+1)/ul(im);
ul=ratio*ul;
u(1:im)=ul(1:im);
y(1:im)=u(1:im).^2;

% Find the wavefunction on the right

u(im+1:nx)=fliplr(ur(1:nx-im));
y(im+1:nx)=u(im+1:nx).^2;

%Normalize the wavefunction
sum = simpson(y, h);
sum = sqrt(sum);
u=u/sum;


% Method to achieve the evenly spaced Simpson rule.
function xx=simpson(y,h)
n = length(y)-1;
s0 = 0; s1 = 0; s2 = 0;
for i=2:2:n
s0=s0+ y(i);
s1 = s1+y(i-1);
s2=s2+ y(i+1);
end
s = (s1+4*s0+s2)/3;
%Add the last slice separately for an even n+1
if mod((n+1),2)== 0% mod
xx= h*(s+(5*y(n)+8*y(n-1)-y(n-2))/12);
else
xx= h*s;
end
   
     
% Method to perform the Numerov integration.
function u=numerov( m,  h, u0,  u1,  q,  s) 
u = ones(1,m);
u(1) = u0;
u(2) = u1;
g = h*h/12;
for i=2:m-1
    c0 = 1+g*q(i-1);
    c1 = 2-10*g*q(i);
    c2 = 1+g*q(i+1);
    d = g*(s(i+1)+s(i-1)+10*s(i));
    u(i+1) = (c1*u(i)-c0*u(i-1)+d)/c2;
end

    
     

