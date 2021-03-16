clear, clc, clf
close all
full_white_value=64; % äldre matlabversion - detta värde plottas som vitt (max) med image-kommandot
%full_white_value=255; % nyare matlabversion - prova denna om din plot verkar mörk!
load T_DOE_gen2
%DOE = load ('T_DOE_gen2.mat')
%DOEmat = cell2mat(struct2cell(DOE));
N=1024; % NxN är matrisstorleken (rekommenderad storlek N=1024)
sidlaengd_Plan1=4e-3; % det samplade områdets storlek (i x- eller y-led) i Plan 1 (rekommenderad storlek 4 mm)
a=sidlaengd_Plan1/N; % samplingsavstånd i Plan 1 (och Plan 2 eftersom vi använder PAS)
L=100e-3; % propagationssträcka (dvs avstånd mellan Plan 1 och 2)

lambda_noll=633e-9; % vakuumvåglängd för rött ljus från en HeNe-laser
n_medium=1; % brytningsindex för medium mellan Plan 1 och 2
k=2 * pi /(lambda_noll) ; % k-vektorns längd *** 

xvekt=-N/2*a:a:(N/2-1)*a; % vektor med sampelpositioner i x-led
yvekt=xvekt; % och y-led
[xmat,ymat]=meshgrid(xvekt,yvekt); % koordinatmatriser med x- och y-värdet i varje sampelposition
rmat=sqrt(xmat.^2+ymat.^2); % avståndet till origo i varje sampelpunkt. Observera att alla operationer är elementvisa!

%******* Fält i Plan 1
f_lins=100e-3; % fokallängd på linsen före Plan 1
T_lins=exp(-1i*k*rmat.^2/(2*f_lins)); % Transmissionsfunktion för en lins (linsen är TOK)

D_apertur=2e-3;
T_apertur=rmat<(D_apertur/2); % Transmissionsfunktion för en cirkulär apertur ("pupill") 

omega_in=1e-3; % 1/e2-radie (för intensiteten, dvs 1/e-radie för amplituden) för infallande Gaussiskt fält
E_in_gauss=exp(-rmat.^2/omega_in^2); % Infallande fält: Gaussiskt med plana vågfronter och normalinfall (dvs konstant fas, här=0)

E_in_konstant=ones(N,N); % Infallande fält: Plan våg med normalt infall

E_in_hermitegauss  =  E_in_gauss.*xmat;

E1_gauss=E_in_gauss.*T_lins; % Fältet i Plan 1 (precis efter linsen) för gaussisk stråle 
E1_cirkular=E_in_konstant.* T_apertur .* T_lins; % Fältet i Plan 1 (precis efter linsen) för konstant fält som passerat genom cirkulär apertur *** Ej klar 
E1_hermitegauss = E_in_hermitegauss.* T_lins;

E1=E1_gauss; % Välj fall!

I1=abs(E1).^2; % intensiteten är prop mot kvadraten på fältets amplitud (normalt struntar man i proportionalitetskonstanten)

C = 1;
Dspot = C * lambda_noll/2* omega_in * L

%Dspot vid Gaussisk stråle
% Dspot mäts upp tlll ca 0.06mm f = L = 10cm
NewG10_G = (6e-5 * 2 * omega_in) / (lambda_noll * L) % C = 1.90

% Dspot mäts upp til ca 0.4mm F=L=1m
NewG100_G = (0.64e-3 * 2 * omega_in) / (lambda_noll * L) % C = 1.26



%Dspot vid cirkulär ståle
% Dspot mäts upp tlll ca 0.08mm f = L = 10cm
NewG10_C = (0.07e-3 * 2 * omega_in) / (lambda_noll * L) % C = 2.53

% Dspot mäts upp til ca 0.8mm F=L=1m
NewG100_C = (0.7e-3 * 2 * omega_in) / (lambda_noll * L) % C = 1.53




figure(1)
image(xvekt*1e3,yvekt*1e3,I1/max(max(I1))*full_white_value)
title(['Intensitet i Plan 1. Verkar OK, eller?'])
xlabel('x[mm]')
ylabel('y[mm]')
colormap(gray) 
drawnow
axis('equal')

figure(2)
imagesc(xvekt*1e3,yvekt*1e3,angle(E1))
title(['Fas i Plan 1. Verkar OK, eller?'])
xlabel('x[mm]')
ylabel('y[mm]')
colormap(gray) 
colorbar
drawnow
axis('equal')


%pause % tryck på valfri tangent för att fortsätta

%**** Och nu propagerar vi till Plan 2!
E2=PAS(E1,L,N,a,lambda_noll,n_medium); % Propagation med PAS-funktionen *** Ej klar
%E2 = E2 + 0.5;
I2=abs(E2).^2; 

mattnadsfaktor_plot=500; % anger hur många gånger maxvärdet ska vara mättat i plotten (>1, kan vara bra om man vill se svagare detaljer)
figure(3)
image(xvekt*1e3,yvekt*1e3,I2/max(max(I2))*full_white_value*mattnadsfaktor_plot)
title(['Intensitet efter ' num2str(L*1e3) ' mm propagation med cirkulärt infallande ljus'])
xlabel('x[mm]')
ylabel('y[mm]')
colormap(gray) 
drawnow
%axis(-1.5,1.5,-1.5,1.5)
txt = 'D_{spot} \approx 0.07 mm';
text(0,0,txt,'FontSize',9)


figure(4)
plot(xvekt*1e3,I2(N/2+1,:))
title(['Intensitet längs x-axeln efter ' num2str(L*1e3) ' mm propagation med cirkulärt infallande ljus.'])
xlabel('x[mm]')
drawnow

L = 20e-3;
f_Vlins=150e-3;
T_Viseslins=exp(-1i*k*rmat.^2/(2*f_Vlins));

E1 = ones(N,N);

f_lins=20e-3; % fokallängd på linsen före Plan 1
T_lins=exp(-1i*k*rmat.^2/(2*f_lins));

E2 = E1 .*T_lins .* T_DOE_gen2 .* T_Viseslins;

E2=PAS(E2,20e-3,N,a,lambda_noll,n_medium);

I2=abs(E2).^2; 
mattnadsfaktor_plot=150; % anger hur många gånger maxvärdet ska vara mättat i plotten (>1, kan vara bra om man vill se svagare detaljer)
figure(5)
image(xvekt*1e3,yvekt*1e3,I2/max(max(I2))*full_white_value*mattnadsfaktor_plot)
camroll(180)
%title(['Ofarliga meddelandet, mattnadsfaktor=' num2str(mattnadsfaktor_plot), 'VerticalAlignment', 'bottom'])

t=title('Title1')
   xlabel('farliga meddelandet, mattnadsfaktor=150, camroll 180\circ')
   set(t,'position',get(t,'position')-[0 1.4 0])
%xlabel('x[mm]')
%ylabel('y[mm]')
colormap(gray) 

drawnow
axis('equal')


%Ofarliga meddelandet "Make Newton great again"
%Farligt meddelande "<3 Fresnel"

