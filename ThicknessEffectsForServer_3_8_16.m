function out = ThicknessEffectsForServer()

hold on
hold all

ratio = 0.9;
density = 980; %kg/m^3
maxratio = ratio;
R = 0.05;%m
mass = (4/3)*pi*R^3*density;%kg
viscosity = 0.001; 
v=-0.003673469;

% drive velocity, m/s
h0=1.08994E-05;%m
time=0;
poi=.5;
E=1e6;% pa
Rs = [0:R/1000:R/10];
YvsR = zeros(200,length(Rs));
YvsR(1,:) = Rs.^2/(2*R) + h0; 
hinitial = YvsR(1,:) ;
tincre = 0.0002;

% Data initialization finished

hpass1=hinitial;
allforce=zeros(400,1);
allhvst=zeros(400,length(Rs));
allwvst=allhvst;
alltvst=zeros(400,1);
allpvst=allhvst;
allevst=zeros(400,length(Rs));
allpevst=zeros(400,length(Rs));
hguess1=zeros(1,length(Rs));% this variable is the guess of new position 
p1=zeros(1,length(Rs));
p1guess=zeros(1,length(Rs));
hcalc=p1guess;
wguess=p1guess;

allhcalc = zeros(3000,length(Rs));
allhguess1 = zeros(3000,length(Rs));

indexer=2;
allhvst(1,:) = YvsR(1,:);
time=time+tincre;
hguess1=YvsR(1,:);
hcalc=YvsR(1,:);
hpreviousSol=YvsR(1,:);
calcint = zeros(1,length(Rs));
dhdt1 = zeros(1,length(Rs)); 
fittinginside = zeros(1,length(Rs));
inside = zeros(1,length(Rs));
ecalc = zeros(1,length(Rs));
peguess = zeros(1,length(Rs));
eguess = zeros(1,length(Rs));
vguess = v;
vpreviousSol = v;
wpreviousSol = zeros(1,length(Rs));

%indexer counts the number of time increments
while indexer<1000 %The outer loop, when the loop execute once one certain time separation was determined

epreviousSol=ecalc; 
hpreviousSol=hcalc;
vpreviousSol=vguess;
wpreviousSol=wguess;
if min(hcalc) < 10e-9
    break
end

criteriacount = 0;
if indexer > 2 % This if statement is for the updating initialize dh.
    hguess1=2*hcalc - allhvst(indexer-2,:);
else
    hguess1 = hinitial + v*tincre;  
%else 
    %hguess1=2*hcalc - allhvst(indexer-1,:);
end 

criteria=false;
runs=1;

while criteria==false% The inner loop
%Solving for dp/dr in lubrication eqn
dhdt1=(hguess1-hpreviousSol)./tincre; 
inside = dhdt1.*Rs;
fittinginside = spline(Rs,inside);
integ = fnint(fittinginside); %indefinite integral of r*dh/dt with respect to r
calcint = fnval(integ,Rs);

for a = 1:length(p1guess)
   dp1finder(a) = -12*viscosity*calcint(a)/(Rs(a)*hguess1(a)^3);
   dp1finder(1) = 0;
end  %dp/dr has been solved

%
for a=length(p1guess)-1:-1:1
   
    p1guess(a)=p1guess(a+1)+((dp1finder(a)+dp1finder(a+1))/2)*(Rs(a+1)-Rs(a));%edge is zero, matrix of dp converted to pressure profile
    eguess(a) = (p1guess(a)+viscosity*epreviousSol(a)/tincre)/(E+viscosity/tincre);
    peguess(a) = eguess(a)*E;
    if peguess(a) > p1guess(a)
        peguess(a) = p1guess(a);
    end

end




p1guess(1) = p1guess(2);
peguess(1) = peguess(2);
I0=pi/2*trapz(Rs,p1guess);
force=trapz(Rs,2*pi*Rs.*p1guess);

for a=1:length(Rs)
%Itail=-p1guess(a)*(pi*Rs(a)^2./(8*Y));%?
dummyp=zeros(1,length(Rs));
dummyp=p1guess;
c=4*Rs(a)*Rs;
d=((Rs+Rs(a)).^2);
phi1=Rs./(Rs+Rs(a)).*ellipke(c./d);
phi1(a)=pi/2;
Ir=trapz(Rs,(dummyp-peguess(a)).*(phi1-pi/2));%?
wguess(a)=4*(1-poi^2)/(pi*E)*(I0+Ir);%

end

% from dummy p =... to wguess(a): determine deformation from pe.
% to update for layering effects, uncomment following code and comment
% above code.

p = peguess;

%for i = 1:inc:maxlength
%Z(1,(i+inc-1)/inc) = -1 * i*trapz(r,r.*p.*besselj(0,i*r));
%end
%A = @(ksi) [(lambda+miu)*ksi miu -1*(lambda+miu)*ksi miu;
    %(lambda+miu)*ksi -1*lambda (lambda+miu)*ksi lambda;
    %-1*ksi.*exp(-1*ksi*h) exp(-1*ksi*h)-h*ksi.*exp(-1*ksi*h) ksi.*exp(ksi*h) exp(ksi*h)+h*ksi.*exp(ksi*h);
    %((lambda+miu)*ksi.*exp(-1*ksi*h))/miu (2+(lambda+miu)*h*ksi/miu).*exp(-1*ksi*h) ((lambda+miu)*ksi/miu).*exp(ksi*h) (-2+(lambda+miu)*h*ksi/miu).*exp(ksi*h)];
%B = @(ksi) [Z(1, (ksi+99)/100)/(2*(ksi.^3)); 0; 0; 0];
%coeff = zeros(4, maxlength/inc);
%format shortEng
%for i = 1:inc:maxlength
%coeff(:, (i+inc-1)/inc) = (A(i))\B(i);
%end
%ratio = zeros(1,length(r));

%for j = 1:length(r)
  %  u = range.*(coeff(1,:).*(range.^2)+coeff(3,:).*(range.^2)-2*coeff(2,:).*range+2*coeff(4,:).*range-(lambda+2*miu)*(coeff(1,:)+coeff(3,:)).*(range.^2)/miu).*besselj(0,range.*r(1,j));
    %wguess(1,j) = trapz(range,u);
   % uLast10 = u(0.9*maxlength/inc:maxlength/inc);
    %ratio(1,j) = max(uLast10)/max(final(1,j));
   % if ratio(1,j) > 0.01
 %plot(range,u)
        %break
    %end
%end

%if indexer==2
    %v = 0.286/0.5*v
%else
    %v = -350E-9
%end

vguess = force*tincre/mass + vpreviousSol;
hcalc = hpreviousSol  + vguess *tincre + wguess-wpreviousSol;

if max(abs(hcalc-hguess1))< 0.5E-10
    criteria=true;
    allruns(indexer,1) = runs;
    ratio = maxratio;
else

allhcalc(runs,:) = hcalc;   
allhguess1(runs,:) = hguess1;
hguess1=hguess1*ratio+hcalc*(1-ratio);   
if indexer > 10
if runs > 5
if (min(hguess1-allhguess1(runs-1)) > 0) && (min(allhguess1(runs-1,:)-allhguess1(runs,:)) > 0)
    ratio = ratio + (1-ratio)*0.1;
    hguess1=allhguess1(runs,:)*ratio+allhcalc(runs,:)*(1-ratio);
    criteriacount = 0;
elseif (max(hguess1-allhguess1(runs-1,:)) < 0) && (max(allhguess1(runs-1,:)-allhguess1(runs,:)) < 0)
    ratio = ratio + (1-ratio)*0.1;
    hguess1=allhguess1(runs,:)*ratio+allhcalc(runs,:)*(1-ratio); 
    criteriacount = 0;
else
    criteriacount = criteriacount +1;
end
maxratio = max(ratio, maxratio);

if (criteriacount > 50) && (max(abs(hcalc-hguess1)) < 3e-10) 
    ratio = (ratio-0.1)/0.9;
    criteriacount = 0;
end
else 
end
end
%if indexer > 50
    % = 0.95;
%end

%if indexer > 57
   % ratio = 0.97;
%end
    

%if indexer > 62
    %ratio = 0.99;
%end

%if indexer > 75
   % ratio = 0.995;
%end

%if indexer > 90
    %ratio = 0.999;
%end

  %  elseif indexer>30

   %     hguess1=hguess1*0.9999+hcalc*0.0001; 

    %elseif indexer>40

       % hguess1=hguess1*0.99999+hcalc*0.00001; 

   % elseif indexer>80

       % hguess1=hguess1*0.999999+hcalc*0.000001;  

   % elseif indexer<20

     %  hguess1=hguess1*0.9999999+hcalc*0.0000001;    
    if rem(runs,10)<1
     subplot(1,4,3)
   plot(Rs,allhcalc(runs,:)-allhguess1(runs,:))
   title(num2str(ratio));
    axis([0 0.001 -1e-9 1e-9])
    drawnow    
      subplot(1,4,2)
   hold on
    plot(Rs,hguess1,'red')
    plot(Rs,hcalc,'blue')
   axis([0 100e-6 0e-9 h0])
    title(num2str(runs))
    drawnow

    end

    runs=runs+1;

    
if runs > 10000
        break
end
    %dhdtfitreal = fitdhdt(1)*Rs.^5+fitdhdt(2)*Rs.^4+fitdhdt(3)*Rs.^3+fitdhdt(4)*Rs.^2+fitdhdt(5)*Rs+ fitdhdt(6);

    %subplot(1,4,4)

    

    %plot(Rs,dhdt1,'*-')

    %hold on

    %fnplt(fittingdhdt1,'red')

    %axis([0 Rs(40) -72E-9 5E-9])

    %hold on

 % drawnow
end

end
cla(subplot(1,4,2))
cla(subplot(1,4,3))
subplot(1,4,1)
hold on
plot(Rs,hcalc)
axis([0 250e-6 0 h0])
title(num2str(indexer));
drawnow




%subplot(1,4,4) 

%hold on 

%plot(Rs,peguess./p1guess)

%axis([0 250e-6 0 1])

%title(num2str(thicknessmax(qw)));

%drawnow 

ecalc=eguess ; 
allforce(indexer,1)=force;
allhvst(indexer,:)=hcalc;
allpvst(indexer,:)=p1guess;
allwvst(indexer,:)=wguess;
allevst(indexer,:)=eguess;
allpevst(indexer,:)=peguess;
allvvst(indexer,1)=vguess;

time=time+tincre;

%plot(Rs,hcalc,'x-')
%axis([0 Rs(20) 0 hinitial*3])
%drawnow
%allhvst(indexer,:)=hcalc;
indexer=indexer+1;

%time=time+tincre;

if runs > 10000
        break
end
end

Stokes = mass * v/(6*pi*viscosity*R^2);
qws = num2str(Stokes);
qwsstring = strcat(qws,'.xlsx');

xlswrite(qwsstring,allhvst,1)
xlswrite(qwsstring,allpvst,2)
xlswrite(qwsstring,allwvst,3)
xlswrite(qwsstring,allforce,4)
xlswrite(qwsstring,allruns,5)
xlswrite(qwsstring,tincre,6)
xlswrite(qwsstring,allvvst,7)

