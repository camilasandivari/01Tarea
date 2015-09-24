import numpy as np
from pylab import *
import astropy as astro
from astropy import units as u
from astropy import constants as const
import scipy as sc
import scipy.integrate
global x
global y
global xmin
global xmax



a=np.loadtxt('sun_AM0.dat')
print type (a)
x=range(len(a))
y=range(len(a))

for i in range(len(a)):
    x[i]=a[i][0]      #guardando l.o
    y[i]=a[i][1]      #guardando espectro

plot(x ,y)
xscale('log')
#yscale('log')
xlabel('longitud de onda [$nm$]')
ylabel('Flujo [$W*m^2*nm^1$]')
title('flujo v/s longitud de onda')
savefig("1.png")
show()



xmin=x[0]
xmax=x[len(x)-1]

def integral(y,x):     #integral para los datos
    integr=0
    i=0
    while i<(len (x)-1):
        paso=x[i+1]-x[i]
        integr += (y[i]+y[i+1])* (paso/2)
        i += 1
    return integr


comparacion1=scipy.integrate.trapz(y,x)
soldatos=integral(y,x)
s= 2.81*(10**23)
lumin= soldatos*s
print ('integracion por trapecio '+str(soldatos))
print ('integracion con funcion de scipy '+str(comparacion1))
print ('luminosidad '+ str(lumin))


h=const.h
c=const.c
k=const.k_B
T=5778  #K

solanalitica=(2*math.pi*h/(c**2))*((k*T/h)**4)*(((math.pi)**4)/15)
print ('solucion analitica planck '+str( solanalitica))

def f(r):           #funcion para planck usando cambio de variable a la tangente
    return (np.tan(r) ** 3 + np.tan(r)**5)/(np.expm1(np.tan(r)))


xmin=0.0000001 #para que no se indefina
xmax= math.pi/2 #por cambio de variable
div=1000

def integral1(f,xmin,xmax,div):  #integral para la funcion
    delta1=np.linspace(xmin,xmax,num=div)
    sep= (delta1[1]-delta1[0])/2
    ii=0
    integrando=0
    while ii < div-2:
        integrando += (sep/3.0)*(f(delta1[ii])+(4*f(delta1[ii+1]))+ f(delta1[ii+2]))
        ii+= 1
    return integrando


P=(2*math.pi*h/(c**2))*((k*T/h)**4)*(integral1(f,xmin,xmax,div))
comparacion2=scipy.integrate.quad(f,xmin,xmax)
P1=(2*math.pi*h/(c**2))*((k*T/h)**4)*comparacion2
a=const.au
b=soldatos/P
c=b**0.5
radio= c*a

print('integracion con funcion scipy'+str(P1))
print ('solucion planck con integracion simson '+ str(P))
print ('radio solar '+str(radio))
#print 'hola'
