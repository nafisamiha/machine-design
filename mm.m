clc;
%specifications
Q=input('Enter the value of kVA rating=\n');%in my case-25+092=117kVA;
VLp=input('Enter the voltage of primary winding=\n');
VLs=input('Enter the voltage of secondary winding=\n');
f=input('Enter frequency=');%may be use to design distribution transformer
Bmax=input('Enter the value of magnetic fluxdensity for CRGO in Wb/m^2=\n');
n=input('Enter the number of steps of core=\n');
ki=input('Enter the value of iron space facor as steps =\n');
ks=input('Enter the stacking factor for paper insulation=\n');
kw=input('Enter Windows space factor=\n');
delta=input('Enter current density in A/mm^2=\n');%forced cooling
r=input('Enter the height width ratio of window Hw/Ww=\n');
Nl=input('Enter the number of layer of lv winding=\n');
dc=input('thickness of rectangular conductor in mm=\n');
wc=input('width of rectangular conductor in mm =\n');
N=input('Enter the number of layer for hv winding=\n');
dsl=input('Enter the distance between core & lv coil in mm=\n');
ds=input('Enter the distance between lv & hv winding in mm=\n');
Nc=input('Enter the number of coil of hv winding=\n');
dbc=input('Enter the distance betweent two coil of hv winding in mm=\n');
p20=input('Enter the resistivity at 20degree C in ohm/mm^2/m=\n');
c20=input('Enter the temperature coeeficient of resistivity at 20 =\n');
Tc=input('Enter the temperature on full load in celsius=\n');
% finding voltage per turn 
k=sqrt(1000/3)/40%no of legs 3 phase core type transformer; 
Et=round(k*sqrt(Q))
%%specific magnetic loading
%Bmax=input('Enter the value of magnetic flux density for CRGO=');
%% cross section area of the core
Ai=(Et*1e6)/(4.44*Bmax*f)%net area of iron core in mm^2
%%diameter of circumference cirle of the core

d=round(sqrt((4*Ai)/(pi*ki*ks)),2,'significant')%diameter
fprintf('checking Bmax');
Ai=ki*ks*(pi*d^2)/4
Bmax=round((Et/(4.44*50*Ai*1e-6)),2,'significant')
%% windows area Aw
%kw=input('Enter Windows space factor=');
%delta=input('Enter current density in A/mm^2=');%forced cooling
Aw=(Q*1e6*1e3)/(3.33*Ai*kw*delta*Bmax*f)
%r=input('Enter the height width ratio of window Hw/Ww=');
t=d:-10:d-100;
for i=length(t)
    Ww=round(((Aw/r)^0.5),2,'significant');
    Wh=round((Aw/Ww),2,'significant'); 
    if (Wh>200)
        Wh=2.*Ww;
         Aw=round(Wh*Ww);
    fprintf('windows width=%.2f\n',Ww);
    fprintf('windows height=%.2f\n',Wh);
     fprintf('windows area=%.2f',Aw);
   end
end
%% main dimension of the core
%diameter=d
Wwc=.95*d%largest width of the core
D=Ww+Wwc%distance between two adjacent limb
Hw=round(Wh+26*2)%height of window
Wwt=round((D*2+Wwc),2,'significant')%total width
Hwt=round((Hw+Wwc+Wwc),2,'significant')%total height
%% noumber of turn in LV & HV winding
 Vs=VLs/sqrt(3);%voltage per phase
 T2=round(Vs/Et)%turn of lv winding
 T11=round(VLp/Et)%turn of hv winding
 t1=-.025:.025:.05%tapping
 for j=length(t1)
         T=T11+t1*T11
        
 end
  T1max=T11+T11*.05
% area of hv & lv winding
I2=(Q*1000)/(sqrt(3)*VLs)%current per phase in lv winding
% helical cylindrical coil is choosen
a2=I2/delta%area of lv conductor
%4.5mm thickness X 6.4mm width and  2 conductor strips are choosen
a2=round((dc*wc)*2,2,'significant')
I1=(Q*1000)/(3*VLp)%current in hv winding & disc coil is choosen
a1=I1/delta%area of hv winding
d1=round((sqrt((a1*4)/pi)),2,'significant')% diameter of conductor in mm
a1=(pi*d1^2)/4
Ac=2*(a1*T1max+a2*T2)%copper area in window
Kwc=Ac/Aw
if (Kwc>.23)
    fprintf('windows space factor is near about=%.2f\n',kw);
end
%% Design and layout of lv winding
fprintf('\nfor paper insultion the size of conductor will be as \n');
dc1=dc+.25;
fprintf('\nthickness=%.2f\n',dc1);
wc1=wc+.25;
fprintf('\nwidth=%.2f\n',wc1);
t2=T2/Nl;%turn per layer
fprintf('\nTurn per layer=%.2f\n',t2);%two layer of conductor is choosen
wdc=(wc1+wc1)*Nl;
fprintf('\nwidth of conductor for two layer=%.2f\n',wdc);
fprintf('\nheight of conductor for each layer=%.2f\n',dc1);
Hwl=round(t2*dc1);
fprintf('\nheight of lv winding in window=%.2f\n',Hwl);
fprintf('\nthickness of lv coil=%.2f\n',wdc);
fprintf('\ndistance between core & lv coil=%.2f\n',ds);
di=d+ds*2;
fprintf('\ninside diameter of lv coil=%.2f\n',di);
do=di+2*wdc;
fprintf('\noutside diameter of lv winding=%.2f\n',do);
dm=di+wdc;
fprintf('\nmean diameter of lv winding=%.2f\n',dm);
lm=pi*dm;
fprintf('\nmean length of turn of lv coil=%.2f\n',lm);
%% Design and layout of hv winding
dih=do+ds*2;
fprintf('\ninside diameter of lv coil=%.2f\n',dih);
Ts=round(T1max/Nc);
fprintf('\nsplit hv winding in 4 coil each with turn=%.2f\n',Ts);
dsh=d1+.25;
fprintf('\nsize of cnductor=%.2f\n',dsh)
Tsn=Ts/N;%turn per layer
fprintf('\nTurn per layer=%.2f\n',Tsn);
Hwh=Tsn*dsh;
fprintf('\nheight of conductor for each hv coil=%.2f\n',Hwh);
wt=N*dsh;
fprintf('\nthickness of each coil=%.2f\n',wt);
doh=dih+2*wt;
fprintf('\noutside diameter of hv coil=%.2f\n',doh);
dmh=dih+wt;
fprintf('\nmean diameter of hv winding=%.2f\n',dmh);
lmh=pi*dmh;
fprintf('\nmean length of turn of hv coil=%.2f\n',lmh);
hhv=Hwh*Nc+3*dbc;
HW=hhv+22*2;
fprintf('\nHeight of Window=%.2f\n',HW);
%% percentage reactance
Lmt=(lm+lmh)/(2*1000);
fprintf('\naverage mean length of turn in meter=%.2f\n',Lmt)
AT=I2*T2;
Hc=(HW+Hwl)/(2*1000);
fprintf('\nmean height of coils in meter=%.2f\n',Hc);
a=ds/1000;%distance between lv & hv;
b1=wt/1000;%width of hv coil;
b2=wdc/1000;%thickness of lv coil;
u0=(4*pi)*10^-7;
 pX=((((2*pi*f*u0*Lmt*AT)/(Hc*Et))*(a+(b1+b2)/3))*100);
 fprintf('% Reactance X=%.2f',pX);
 %% percentage resistance
 %resistance of lv winding:(per phase)
 p75=p20*(1+c20*(Tc-20));
 Rlv=(p75*lm*T2)/(a2*1000);
 fprintf('resistance of lv winding=%.5f\n',Rlv);
 %resistance of hv winding:(per phase)
 Rhv=(p75*lmh*T11)/(a1*1000);
 fprintf('resistance of hv winding=%.5f\n',Rhv);
 %equivalent resistance referred to h.v winding(per phase)
 TR=(VLp*sqrt(3))/200;
 R=Rhv+Rlv*TR^2;
 fprintf('equivalent resistance referred to h.v winding=%f\n',R)
 pR=(((I1*R)/VLp)*100);
 fprintf('% Resistance=%.2f\n',pR)
 %% percentage impedance
 pZ=sqrt(pX.^2+pR.^2);
 fprintf('PERCENTAGE IMPEDENCE=%.2f',pZ);
 
 
 
 
 
 
 








  



 
 
 
 
 
 













