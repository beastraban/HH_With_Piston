function [T,N,M,H,VV,PP,XU,XD,Fext]=HHmodelWpiston(gkFactor,gmFactor,ghFactor)
%Best results in: 
%50,5,0.5 which come out to be 0.103, 0.0103, 0.00103

Cm=1;
%R=
an=@(x) 0.01*(x+55.0)/(1.0 - exp(-(x+55.0) / 10.0));
am= @(x) 0.1*(x+40.0)/(1.0 - exp(-(x+40.0) / 10.0));
ah=@(x) 0.07*exp(-(x+65.0) / 20.0);

bn=@(x) 0.125*exp(-(x+65) / 80.0);
bm=@(x) 4.0*exp(-(x+65.0) / 18.0);
bh=@(x) 1.0/(1.0 + exp(-(x+35.0) / 10.0));


Vk=-77;
Vm=50;
Vl=-54.387;

Dk=1;
Dm=1;
Dl=1;

gk= 36;%*(1-Dk*x);    %K
gm= 120;%*(1-Dm*x);   %Na  
gl= 0.3;%*(1-Dl*x);   %Leak

Iext=@(x) 0;
Vin=-65;
for i=1:20
    Iin=@(x) 7*(heaviside(x-(2*i-1)*500/40)-heaviside(x-(2*i)*500/40));%+30*(heaviside(x-300)-heaviside(x-400)
    Iext=@(x) Iext(x)+Iin(x);
end
%define Iext
%*(heaviside(x-100)-heaviside(x-120))%+30*(heaviside(x-300)-heaviside(x-400));

%define external force
D=@(x) 0.04*(0.5*((tanh((x-100)/10))-(tanh((x-180)/10))));
tt=0:0.01:500;
figure()
plot(tt,Iext(tt));

Vdout=zeros(size(tt));
FF=zeros(size(tt));
P0=30;
L=1;
A=1;
Rin=0;%0.1;
Rout=1;
Ku=1;
Kd=1;
Mu=1;
Md=1;
Pout=30.001;
Pin=29.999;

y0=[0.32,0.05,0.6,Vin,P0,0,0,0,0,]; %n,m,h,V,P0,Xup0,Xud,Xdown0,Xdd
y(:,1)=double(y0');
t=tt(1);
yy=y;
T=tt.';
opts=odeset('OutputFcn',@odeplot);


counter=1;    
factor=3.5;  
[t,yy]=ode45(@fm,T,yy);%,opts);





function YOUT=fm(t,Y)
        n=Y(1);
        m=Y(2);
        h=Y(3);
        V=Y(4);
        P=Y(5);
        Xu=Y(6);
        Xud=Y(7);
        Xd=Y(8); 
        Xdd=Y(9);
        
       
        dndt=an(V)*(1-n)-bn(V)*n;
        dmdt=am(V)*(1-m)-bm(V)*m;
        dhdt=ah(V)*(1-h)-bh(V)*h;
        %Iext(t)
        Xddd=(A*(P-P0)-Kd*Xd -Xdd)/Md;
        if (Xd>0.003)&&(Xddd>0)
            Xddd=0;
            Xd=0.003;
            Xdd=0;
        end
        
        
        
        Xudd=-((A*(P-P0))+Ku*(Xu-D(t)) +Xud)/Mu;
        Vol=A*(L-Xu+Xd);
        Pd=-(P/Vol)*A*(Xdd-Xud)-Rout*heaviside(P-Pout)*(P-Pout)+Rin*heaviside(Pin-P)*(Pin-P);
        
        if (Xu>0.03)&&(Xud>0)
            %Xudd=0;
            Xu=0.03;
            Xud=0;
        end
        
        Vd= (-gk*(1-Xd*gkFactor)*n^4*(V-Vk) - gm*(1-Xd*gmFactor)*m^3*h*(V-Vm) - gl*(1-Xd*ghFactor)*(V-Vl) +Iext(t))/Cm;
        Vdout(counter)=Vd;
        FF(counter)=D(t);
        counter=counter+1;
        %Vdd=I+Iext;
        %We set R=1 and a=1 for simplicity.
        YOUT=[dndt,dmdt,dhdt,Vd,Pd,Xud,Xudd,Xdd,Xddd].';
end


T=t;    
N=yy(:,1);
M=yy(:,2);
H=yy(:,3);
VV=yy(:,4);
VVd=Vdout.';
PP=yy(:,5);
XU=yy(:,6);
XD=yy(:,8);
Fext=FF.';
figure()
subplot(2,1,1)
plot(T,VV);
hold on
subplot(2,1,2)
plot(T,XD)
hold on
max(XD)
max(XD*gkFactor)
max(XD*gmFactor)
max(XD*ghFactor)
figure()
plot(T,FF)
end
