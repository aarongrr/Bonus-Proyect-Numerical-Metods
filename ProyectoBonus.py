""" Proyecto bonus
Por Javier Aarón Guerrero Alvarez 
Métodos numéricos"""

import math as m
import numpy as np
import sympy as sy

print('Programa para encontrar los ceros de una función ingresada por el usuario.')
print('   ')
x=sy.symbols('x')
f=input('Escribe la expresión: ')
#print(sy.sympify(f).subs(x,m.pi))
df=sy.diff(f,x)
#print(sy.sympify(df))
print('La espresión ingresada es: ')
print(sy.sympify(f))
#TOL=float(input('Ingresar Tolerancia deseada: '))
#M=int(input('Ingresar el número máximo de iteraciones: '))
#--------------------------------------------- 
#DEFINICIÓN DE LOS MÉTODOS
#---------------------------------------------
    
        
opcion = 0
while True:
    print("""
    Menú de métodos disponibles
    
    1) Bisección
    2) Newton
    3) Secante
    4) Newton Modificado
    5) Steffensen
    """)
    opcion = int(input("Elige una opción: ") )     

    if opcion == 1:
        #Bisección
        def s(t):
            return sy.sympify(f).subs(x,t)

        print('Resolviendo s(t) por método de bisección.')
        print('Por favor ingrese el intervalo de búsqueda.')
        a=float(input('a='))
        b=float(input('b='))

        print('s(a)= %6.4f,   s(b)= %6.4f'%(s(a),s(b)))

        while np.sign(s(a)) == np.sign(s(b)):
            print('Por favor ingrese un intervalo que contenga la solución.')
            a=float(input('a='))
            b=float(input('b='))
            print('s(a)= %6.4f,   s(b)= %6.4f'%(s(a),s(b)))

        TOL=float(input('TOL='))
        M=int(input('M='))

        pn=(a+b)/2
        fn=s(pn)

        print('n \t an \t\t bn \t\t pn \t\t f(pn)')
        print('------------------------------------------------------')

        N=1
        while N<M:
            print('%d \t %10.10lf \t %10.10lf \t %10.10lf \t %10.10lf'%(N,a,b,pn,fn))
            if fn==0 or (b-a)/2<TOL:
                break
            if np.sign(s(a))*np.sign(s(pn))>0:
                a=pn
            else:
                b=pn
            pn=(a+b)/2
            fn=s(pn)
            N+=1

        if(N==M):
            print('Fallo. Aumente M o disminuya TOL.')
            
    elif opcion == 2:
        #Newton
        def gx(t):
            return t-((sy.sympify(f).subs(x,t))/(sy.sympify(df).subs(x,t)))

        TOL = float(input('Tolerancia: '))
        M = int(input('Maximo de iteraciones: '))
        p0 = float(input('Punto inicial P0: '))

        #Imprimir tabla
        print('No. \t P1')
        print('---------------')
        print('0 \t {0}'.format(p0))

        #Implementamos el método
        N=int(1)
        while N<M:
            p=gx(p0)
            print('{0} \t {1}'.format(N,p))
            if m.fabs(p-p0)<TOL:
                break
            else:
                p0=p
                N += 1
                
        #Checando si fallo
        if N==M:
            print('Proceso fallido.')
            print('El cero de la funcion no fue encontrado despues de %d iteraciones.'%(M))
            print('Intente aumentar el numero de iteraciones o disminuir TOL.')
        
    elif opcion == 3:
        #Secante
        def fx(t):
            return sy.sympify(f).subs(x,t)

        # Pedir al usuario ingresar p0, TOL y M.
        TOL = float(input('Tolerancia: '))
        M = int(input('Maximo de iteraciones: '))
        p0 = float(input('P0: '))
        p1 = float(input('P1: '))

        f0=fx(p0)
        f1=fx(p1)

        # Imprimir tabla
        print('N \t p1')
        print('---------------')

        print('0 \t {0}'.format(p0))
        print('1 \t {0}'.format(p1))
        # print('0 \t %10.10f'%(p0))

        N=int(2)

        while N<M:
            p=p1-f1*(p1-p0)/(f1-f0)
            print('{0} \t {1}'.format(N,p))
            if m.fabs(p-p1)<TOL:
                break
            N+=1
            p0=p1
            f0=f1
            p1=p
            f1=fx(p)
                
                
        # Checando si fallo
        if N==M:
            print('Proceso fallido.')
            print('El cero de la funcion no fue encontrado despues de %d iteraciones.'%(M))
            print('Intente aumentar el numero de iteraciones o disminuir TOL.')
    
    elif opcion == 4:
        #Newton Modificado
        def gxm(t):
            return t-((sy.sympify(f).subs(x,t))/(sy.sympify(df).subs(x,t)))

        TOL = float(input('Tolerancia: '))
        M = int(input('Maximo de iteraciones: '))
        p0 = float(input('P0: '))

        print('N \t p1')
        print('---------------')
        print('0 \t %1.4e'%(p0))

        N=int(1)
        while N<M:
            p=gxm(p0)
            print('%d \t %1.4e'%(N,p))
            if m.fabs(p-p0)<TOL:
                break
            else:
                p0=p
                N += 1 # N=N+1

        # Checando si fallo
        if N==M:
            print('Proceso fallido.')
    
    elif opcion == 5:
        #Steffensen
        def Stff_x(t):
            return ((sy.sympify(f).subs(x,(t-sy.sympify(f).subs(x,t))))/(sy.sympify(f).subs(x,t)))-1

        print('Metodo de Steffensen')

        p0=float(input('P0 = '))
        TOL=float(input('TOL = '))
        M=int(input('M = '))

        print('N \t Pn')
        print('---------------')
        print('0 \t %10.5f'%p0)

        N=1
        while N<M:
            p1=Stff_x(p0) #Generamos p1 y p2
            p2=Stff_x(p1)
            p=p0-pow(p1-p0,2)/(p2-2*p1+p0) #Aitken para p0, p1, p2
            print('%d \t %10.5f'%(N,p))
            if(m.fabs(p-p0)<TOL):
                break
            else:
                p0=p #Reemplazamos p0
                N=N+1

        if(N==M):
            print('Fallo en el numero max. de iteraciones.')
        sleep(10)
    else:
        break