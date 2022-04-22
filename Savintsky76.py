from math import*
from sympy import*


D=60000              #Desplazamiento en libras
LCG=29               #LCG en [ft]
VCG=2                #VCG en [ft]
B=14                 #B manga [ft]
V=67.5               #Vel en [ft/s]
bet=10               #Angulo de astilla muerta [grados]
fpp=0.5              #Distance [ft] between the propeller axis and the center of gravity, measured normal to the propeller axis
epsilon=0.0698132    #the inclination of the propeller axis relative to the keel
A=1.25               #Facto de interferencia en caso del catamaran
r=0.851              #Facto de separacion cascos en caso de catamaran



Pesp=62.4                                  #Peso espesifico sistema ingles lb/ft3
g=32.2                                     #Gravedad sistema Ingles ft/s2
rho=Pesp/g
visc=(3.281**2)*1.003E-6                   #Viscocidad cinematica ft2/s

tao=[radians(i) for i in range(2,16,1)]    #Angulos de trimado

Cv=V/(sqrt(g*B))                           #Coeficente de Velocidad

Clb=D/(0.5*rho*pow(V*B,2))                 #Coefiente de Fuerza de Sustentacion con bet

#Calculo de Cl0                             Coeficiente de Sustentacion con angulo de astilla muerta =0
x=Symbol('x')
f=x-(0.0065*bet*(x**0.6))-Clb
df=diff(f,x)
r=Clb
for i in range(1,50,1):
    r=r-float(f.subs(x,r))/float(df.subs(x,r))
Cl0=r
    
def DragCal(D,B,V,bet,rho,visc,tao,Cv,Cl0):
    
    A=1.25               #Facto de interferencia en caso del catamaran
    r=0.851              #Facto de separacion cascos en caso de catamaran
    
    x=Symbol('x')                                                                          #Calculo de landa
    f=(( degrees(tao))**1.1)*((0.012*(x**0.5))+((0.0055*(x**(5/2)))/(Cv**2)))-Cl0
    #f=((r**(3/2))*((degrees(tao))**1.1)*((0.012*(x**0.5)/A)+((0.0055*(x**(5/2))*A)/((Cv**2)*r))))-Cl0
    df=diff(f,x)
    r=Cl0
    for i in range(1,50,1):
        r=r-float(f.subs(x,r))/float(df.subs(x,r))
    landa=r
    
    Vm=V*sqrt(1-((0.012*((degrees(tao))**1.1))/(sqrt(landa)*cos(tao))))                  #Velocidad del agua relativa al casco
    
    Re=Vm*landa*B/visc                                                                   #Numero de Reynold
    
    Cf=0.0004+(0.075/pow((log10(Re)-2),2))                                               #Coeficiente de de resistencia friccional con rugocidad de ITTC
    
    Df=rho*((Vm*B)**2)*landa*Cf/(2*(cos(radians(bet))))                                  #Resistencia por friccional [N]
    
    Dp=D*tan(tao)                                                                        #Arrastre componete de presion [N]
    
    Drag=(Df/(cos(tao)))+Dp                                                              #Arrastre total [N]
        
    return (landa,Vm,Re,Cf,Df,Dp,Drag)

vfspp=0 #Contador
vfsnn=0 #Contador
j=0     #Contador


while ((vfsnn+vfspp)<2):
    [landa,Vm,Re,Cf,Df,Dp,Drag]=DragCal(D,B,V,bet,rho,visc,tao[j],Cv,Cl0)                                  #funcion para determinar el arrastre total

    Cp=0.75-(1/((5.21*(Cv**2)/(landa**2))+(2.39)))                                                         #Calculo de Centro de presion

    c=LCG-(Cp*landa*B)                                                                                     #Calculo de c

    a=VCG-(B*tan(radians(bet))/4)                                                                          #Calculo de a

    taof=(D*(((1-(sin(tao[j])*sin(tao[j]+epsilon)))*c/(cos(tao[j])))-(fpp*sin(tao[j]))))+Df*(a-fpp)        #Evaluacion de funcion para evaluar el equilibrio en funcio de tao (trimado)
    #Identificando el cambio de signo de la funcion taof
    if (taof<0):
        vfspp=1
        am=degrees(tao[j])
        fm=taof
    if (taof>0):
        vfsnn=1
        amm=degrees(tao[j])
        fmm=taof
    j=j+1

#Interpolando tao de equilibrio
taoeq=am-(fm*((amm-am)/(fmm-fm)))
print(taoeq," --->Tao equilibrio (Trim)")
taoeq=radians(taoeq)
#Evalua el arrastre total con el tao de equilibrio
[landa,Vm,Re,Cf,Df,Dp,Drag]=DragCal(D,B,V,bet,rho,visc,taoeq,Cv,Cl0)
print(Drag," --->Drag")
#Evalua el EHP
EHP=Drag*V/550
print(EHP," --->EHP")


