clear all;
close all;
clc;
%% Specifications:
  S=79;      %kVA rating
  Vp=6600;    %voltage of primary winding
  Vs=220;     %voltage of secondary winding
  f=50;       %frequency
  n=7;        %number of steps of the core
  Bmax=1.7;   %magnetic flux density for CRGO in Wb/m^2
  delta=2.5;  %current density in A/mm^2%for natural cooling
  r=2;%input('Enter the height width ratio of window Hw/Ww=\n');
  tc=5.1;%input('thickness of rectangular conductor in mm=\n');
  wc=8.1;%input('width of rectangular conductor in mm =\n');
  N1=input('Enter the number of layer for hv winding=\n');
  N2=input('Enter the number of layer for lv winding=\n');
  Nc=input('Enter the number of coil in hv winding=\n');
  dcl=input('Enter the distance between core & lv coil in mm=\n');
  dlv=input('Enter the distance between lv & hv winding in mm=\n');
  d2c=input('Enter the distance between two coil of hv winding in mm=\n');
  %p20=input('Enter the resistivity at 20degree C in ohm/mm^2/m=\n');
  %c20=input('Enter the temperature coeeficient of resistivity at 20 =\n');
  Tc=75; %input('Enter the temperature on full load in celsius=\n');
%% Voltage per turn(Et): 
 %no of legs 3 phase core type transformer 
  Et=((sqrt((S*1000)/3))/40);
  fprintf('voltage per turn,Et=%f\n',Et);

%% Specific magnetic loading:
  fprintf('magnetic flux density for CRGO in Wb/m^2,Bmax=%f\n',Bmax);

%% Cross section area of the core(Ai):
  Ai=((Et*1e6)/(4.44*Bmax*f)); %net area of iron core in mm^2
  fprintf('cross section area,Ai=%f\n',Ai);

%% Diameter of the circumscribing cirle for the core(d):
  ki=.88  %core space factor
  ks=.92  %stacking factor for paper insulation
  x=input('multipling factor=\n');
  d=x*round(sqrt((4*Ai)/(pi*ki*ks)));%,2,'significant')
  fprintf('diameter,d=%f\n',d);
  d1=130;
  fprintf('checking Bmax');
  Ai=ki*ks*(pi*(d1)^2)/4
  Bmax=(Et/(4.44*50*Ai*1e-6))

%% windows area Aw:
  kwi=0.27  %Windows space factor
  Aw=(S*1e6*1e3)/(3.33*Ai*kwi*delta*Bmax*f)
  t=d1:-10:(d1)-100;
  for i=length(t)
      Ww=round((Aw/r)^.5);
      Wh=round(Aw/Ww);
      if(Wh>200)
          Wh=2*Ww;
          Aw=Wh*Ww;
          fprintf('windows width,Ww=%.2f\n',Ww);
          fprintf('windows height,Wh=%.2f\n',Wh);
          fprintf('windows area,Aw=%.2f',Aw);
      end
  end

%% The main dimension of the core:
 %diameter=d1
  Wwc=.78*d1          %largest width of the core
  D=Ww+Wwc           %distance between two adjacent limb
  Hw=(280+17*2); %height of window
  Wtw=round(D*2+Wwc);   %total width
  Hth=round(Hw+Wwc+Wwc);%total height
  fprintf('height of the window,Hw=%.2f\n',Hw);
  fprintf('total width,Wtw=%.2f\n',Wtw);
  fprintf('total height,Hth=%.2f\n',Hth);

%% Number of turns in LV winding;
  Vss=Vs/sqrt(3);    %voltage per phase
  T2=round(Vss/Et+.5)+1;  %turns of lv winding
  fprintf('turns per phase on low voltage winding,T2=%.2f\n',T2);
 
%% Number of turns in HV winding:
  T1=round(Vp/Et);   %turn of hv winding
  fprintf('turns per phase on high voltage winding,T1=%.2f\n',T1);
  t1=-.05:.025:.05   %tapping
  for j=length(t1)
           T=T1+t1*T1;
  end
  T1max=T1*1.05

%% Area of lv winding:
  I2=(S*1000)/(sqrt(3)*Vs) %current per phase in lv winding  
  %helical cylindrical coil is choosen
  a2=I2/delta;             %area of lv conductor
  a2=round((tc*wc)*2);
  fprintf('area of low voltage conductor,a2=%.2f\n',a2);
 
%% Area of hv winding:
  % disc coil is choosen
  I1=(S*1000)/(3*Vp)%current in hv winding per phase
  a1=I1/delta %area of hv winding
  %round conductor is choosen
  D1=round((sqrt((a1*4)/pi)); %diameter of conductor in mm
  D2=1.4;
  a1=(pi*(D2)^2)/4;
  fprintf('area of high voltage conductor,a1=%.2f\n',a1);
  Ac=2*(a1*T1max+a2*T2)    %copper area in window
  Kwf=Ac/39200;              %window space factor
  if (Kwf>.23) 
      fprintf('window space factor is near about,kwi=%.2f\n',kwi);
  end

%% Design and layout of lv winding
  fprintf('\nfor paper insulation the size of conductor will be as \n');
  tc1=tc+.25;
  fprintf('\nthickness,tc1=%.2f\n',tc1);
  wc1=wc+.25;
  fprintf('\nwidth,wc1=%.2f\n',wc1);
  %two layer of conductor is choosen
  t2=T2/N2; %turn per layer
  fprintf('\nTurn per layer,t2=%.2f\n',t2);
  Hlw=round(wc1*t2);
  Tlc=(tc1+tc1)*N2;
  fprintf('\nheight of the lv winding in window,Hlw=%.2f\n',Hlw);
  fprintf('\nthickness of lv coil,Tlc=%.2f\n',Tlc);
  fprintf('\ndistance between core & lv coil,dc1=%.2f\n',dcl);
  Di=d+dcl*2;
  fprintf('\ninside diameter of lv coil,Di=%.2f\n',Di);
  Do=Di+2*Tlc;
  fprintf('\noutside diameter of lv winding,Do=%.2f\n',Do);
  Dm=Di+Tlc;
  fprintf('\nmean diameter of lv winding,Dm=%.2f\n',Dm);
  Lm=pi*Dm;
  fprintf('\nmean length of turn of lv coil,Lm=%.2f\n',Lm);

%% Design and layout of hv winding:
  fprintf('\ndistance between lv & hv coil,d1v=%.2f\n',dlv);
  Dih=Do+dlv*2;
  fprintf('\ninside diameter of hv coil,Dih=%.2f\n',Dih);
  Ts=round(T1max/Nc);
  fprintf('\nsplit hv winding in Nc coil each with turn,Ts=%.2f\n',Ts);
  Dsh=D1+.25;
  fprintf('\nsize of cnductor,Dsh=%.2f\n',Dsh);
  t1h=Ts/N1; %turn per layer
  fprintf('\nTurn per layer,t1h=%.2f\n',t1h);
  Hwh=t1h*Dsh;
  fprintf('\nheight of winding in each hv coil,Hwh=%.2f\n',Hwh);
  Th=N1*Dsh;
  fprintf('\nthickness of each coil,Th=%.2f\n',Th);
  Doh=Dih+2*Th;
  fprintf('\noutside diameter of hv coil,Doh=%.2f\n',Doh);
  Dmh=Dih+Th;
  fprintf('\nmean diameter of hv winding,Dmh=%.2f\n',Dmh);
  Lmh=pi*Dmh;
  fprintf('\nmean length of turn of hv coil,Lmh=%.2f\n',Lmh);
  Hh=Hwh*Nc+3*d2c;
  fprintf('\nHeight of the hv coil in window,Hh=%.2f\n',Hh)
  HW=Hh+26*2;
  fprintf('\nHeight of the  Window,HW=%.2f\n',HW);

%% Percentage reactance:
  Lmt=(Lm+Lmh)/2;
  fprintf('\naverage mean length of turn in mm,Lmt=%.2f\n',Lmt);
  AT=I2*T2;
  Hm=(Hh+Hlw)/2;
  fprintf('\nmean height of coils in mm,Hm=%.2f\n',Hm);
  a=dlv;   %distance between lv & hv;
  b1=Th;   %width of hv coil;
  b2=Tlc;  %thickness of lv coil;
  u0=(4*pi)*10^-7;
  pX=((((2*pi*f*u0*Lmt*AT*10^-3)/(Hm*Et))*(a+(b1+b2)/3))*100);
  fprintf('percentage Reactance,X=%.2f\n',pX);

%% percentage resistance:
 %resistance of lv winding:(per phase)
  p20=0.01724; %ohm/mm^2/m
  c20=0.00393;
  pTc=p20*(1+c20*(Tc-20))
  Rlv=(pTc*Lm*T2)/(a2*1000);
  fprintf('resistance of lv winding,Rlv=%.5f\n',Rlv);
 %resistance of hv winding:(per phase)
  Rhv=(pTc*Lmh*T1)/(a1*1000);
  fprintf('resistance of hv winding,Rhv=%.5f\n',Rhv);
 %equivalent resistance referred to h.v winding(per phase)
  TR=(Vp*sqrt(3))/Vs;
  R=Rhv+Rlv*TR^2;
  fprintf('equivalent resistance referred to h.v winding,R=%f\n',R)
  pR=(((I1*R)/Vp)*100);
  fprintf('percentage Resistance,R=%.2f\n',pR);
%% Percentage impedance:
  pZ=sqrt(pX.^2+pR.^2);
  fprintf('percentage impedence,Z=%.2f\n',pZ);
  
%% Weight of iron in core and yoke assembly:  
  volume=Ai*(Wtw*2+HW*3);
  density_iron=7.85*1000; %in kg/m^3
  Wc_y=(volume*density_iron)/1000^3;
  fprintf('Weight of core and yoke,Wc&y=%.2f\n',Wc_y);

%% Calculation of core or iron loss:
  %Core loss at Bmax=1.70 wb/m^2 is 1.3 watts/kg
  core_loss=(Wc_y)*1.3;
  fprintf('Core loss in transformer,core_loss=%.2f\n',core_loss);
  
%% Magnetising volt amperes:
 % For Bmax=1.70 wb/m^2, VA/kg from the curve is 15 VA/kg
  Mva=Wc_y*15;
  fprintf('Magnetising volt amperes,Mva=%.2f\n',Mva);
 
%% Weight of lv winding:
  %density of copper 8.89 g/cm^3
  %number of turn=T2, area=a2
  %mean lengh of turn=Lm
  Wlv=(8.89*T2*Lm*a2)/1000^2; %per limb
  fprintf('weight of lv winding (per limb),Wlv=%.2f\n',Wlv);

%% Weight of lv winding:
  %number of turns=T1max, normal=T1,area=a1
  %mean lengh of turn=Lmh
  %Weight of 4 coils(one limb):
  W4c=(8.89*a1*Lmh*T1max)/1000^2;
  fprintf('weight of 4 coil(one limb),W4c=%.2f\n',W4c);
  %For normal turns
  W_c=(8.89*a1*Lmh*T1)/1000^2;
  fprintf('weight of the coils(one limb)=%.2f\n',W_c);
  total_weight=3*(Wlv+W_c); %in kg
  fprintf('total weight of copper in transformer,total_weight=%.2f\n',total_weight);
  
%% Copper loss and load loss at 75:
  %hv current per phase= I1
  copper_loss=3*I1^2*R; %copper loss for 3 phases
  fprintf('copper loss for 3 phases,copper_loss=%.2f\n',copper_loss);
  %add stray load loss about 7%
  load_loss=copper_loss*1.07; %at 75%
  total_loss=core_loss+load_loss;
  fprintf('total loss in transformer,total_loss=%.2f\n',total_loss);
  
%% Calculation for performance:
  %efficiency on full load at unity power factor 
  Output=S*1000; %in watts
  efficiency1=(Output*100)/(Output+total_loss);
  fprintf('efficiency on full load at unity power factor,effiency(1)=%.2f\n',efficiency1);
  %efficiency on 3/4th load at unity power factor 
  load_loss(34)=load_loss*(3/4)^2;
  total_loss(34)=core_loss+load_loss(34);
  efficiency2=(.75*Output*100)/(.75*Output+total_loss(34));
  fprintf('efficiency on 3/4th load at unity power factor,effiency(3/4)=%.2f\n',efficiency2);
  %efficiency on 1/2th load at unity power factor 
  load_los=load_loss*(1/2)^2;
  total_los=core_loss+load_los;
  efficiency3=(0.5*Output*100)./(0.5*Output+total_los);
  fprintf('efficiency on 1/2th load at unity power factor,effiency(1/2)=%.2f\n',efficiency3);
  
%% Regulation:
  %Regulation on full load at unity power factor 
  V=1; %in volt
  E=sqrt((V+(pR/100))^2+(pX/100)^2); %in volt
  regulation1=(E-V)*100
  %Regulation on full load at .8 power factor lagging 
  regulation2=(pR*cos(.8)+pX*sin(.6))
 
%% Core loss current,magnetising current,no load current:
  Ic=core_loss/(3*6600);
  fprintf('core loss current,Ic=%.2f\n',Ic);
  Im=Mva/(3*6600);
  fprintf('magnetising current,Im=%.2f\n',Im);
  I0=sqrt(Ic^2+Im^2);
  fprintf('no load current per phase,I0=%.2f\n',I0);
  I0f=I0*100/I1;
  fprintf('no load current w.r.t full load current in parcentage,I0(w.r.t_fl)=%.2f',I0f);
  
  

  
 



  
    
  
  


