# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 18:17:17 2018

@author: Sebastian Cachaya
"""

    
class Derivada:
    def __init__(self, f, metodo , dx= 0.001):
        self.f = f
        self.dx = dx
        self.met=metodo

    def calc(self,x):
        if(self.met == "adelante"):
            return (self.f(x + self.dx) - self.f(x))/self.dx
        elif(self.met == "central"):
            return (self.f(x + (self.dx / 2)) - self.f(x - (self.dx / 2)))/self.dx
        elif(self.met == "extrapolada"):
            f1 = Derivada(self.f, "central", self.dx)
            f2 = Derivada(self.f, "central", (self.dx)/2)
            ff=(4*f2.calc(x) - f1.calc(x))/3
            return ff 
        elif(self.met == "segunda"):
            return (self.f(x + self.dx) + self.f(x - self.dx) - 2*self.f(x))/(self.dx*self.dx)

class Zeros:
    def __init__(self, f, metodo, error=1e-4, max_iter=100):
        self.f = f
        self.error = error
        self.max_iter = max_iter
        self.met=metodo
        
    def zero(self,vi):
        iteracion = 0
        if(type(vi) == type(1.0)):
            error = abs(self.f(vi))
            x0 = vi
            df = Derivada(self.f,"extrapolada",0.00000001)

            if(self.met == "newton"):
                while((error >= self.error) and (iteracion <= self.max_iter)):
                    x0 = x0 - (self.f(x0)/df.calc(x0))
                    error = abs(self.f(x0))
                    iteracion += 1

                return x0
                
            elif(self.met == "newton-sp"):
                root = optimize.newton(self.f,vi)
                return root

            elif(self.met == "fsolve-sp"):
                root = optimize.fsolve(self.f,vi)
                return root

        elif(type(vi) == type((3,14159))):
            if((self.f(vi[0]) < 0) and (self.f(vi[1]) > 0)):
                x1 = vi[0]
                x2 = vi[1]
            elif((self.f(vi[1]) < 0) and (self.f(vi[0]) > 0)):
                x1 = vi[1]
                x2 = vi[0]
                
            if(self.met == "bisectriz"):
                x3 = (x1+x2)/2
                error = abs(self.f(x3))
                iteracion += 1

                while((error >= self.error) and (iteracion <= self.max_iter)):
                    if(self.f(x3) < 0):
                        x1 = x3
                    elif(self.f(x3) > 0):
                        x2 = x3
                    elif(self.f(x3) == 0):
                        return (x3,self.f(x3))
                    x3 = (x1+x2)/2
                    error = abs(self.f(x3))
                    iteracion += 1

                return x3
                
            elif(self.met == "interpolacion"):
                x3 = ((x2*self.f(x1)) - (x1*self.f(x2)))/(self.f(x1)-self.f(x2))
                error = abs(self.f(x3))
                iteracion += 1
                while((error >= self.error) and (iteracion <= self.max_iter)):
                    if(self.f(x3) < 0):
                        x1 = x3
                    elif(self.f(x3) > 0):
                        x2 = x3
                    elif(self.f(x3) == 0):
                        return (x3,self.f(x3))
                    x3 = ((x2*self.f(x1)) - (x1*self.f(x2)))/(self.f(x1)-self.f(x2))
                    error = abs(self.f(x3))
                    iteracion += 1
                return x3

            elif(self.met == "brentq-sp"):
                root = optimize.brentq(self.f,vi[0],vi[1])
                return root

if __name__ == "__main__":
    import math as math
    from scipy import optimize
    f1 = math.sin
    x0=math.pi
    print("\t Metodos para hallar Derivadas:")
    dx0 = Derivada(f1,metodo="adelante",dx=0.001)
    print("Metodo = ",dx0.met,"= ",dx0.calc(x0))
    dx1 = Derivada(f1,metodo="central",dx=0.001)
    print("Metodo = ",dx1.met,"= ",dx1.calc(x0))
    dx2 = Derivada(f1,metodo="extrapolada",dx=0.001)
    print("Metodo = ",dx2.met,"= ",dx2.calc(x0))
    dx3 = Derivada(f1,metodo="segunda",dx=0.001)
    print("Metodo = ",dx3.met,"= ",dx3.calc(x0))
    print("\t Metodos para hallar Raices:")
    x2 = []
    y2 = []
    j = 0
    error = 1e-10
    max_iter = 100
    xi = ((3*math.pi)/4, (5*math.pi)/4)
    while(j<6):
        metodo = "newton" if j==0 else "bisectriz" if j==1 else "interpolacion" if j==2 else "newton-sp" if j==3 else "fsolve-sp" if j==4 else "brentq-sp" 
        raiz = Zeros(f1,metodo,error,max_iter)
        if(j==0 or j==3 or j==4):
            xr = raiz.zero(xi[0])
        else:
            xr = raiz.zero(xi)
        x2.append(xr)
        print("Método = ",metodo,"\t Raíz = ", xr)
        j += 1
