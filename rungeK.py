""" 
McKenna Branting
This is my own work
Runge Kutta Method for solving ODE
"""
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math

t0 = 0
x0 = 1
y0 = 5
h = 0.02
trueSol = 0
Table = np.array([x0,y0,trueSol])
#prints starting values 
print("RUNGE KUTTA METHOD ODE SOLVER")
print("")
print("-------------------------------")

print ("n | x{n} | y{n} | k1 | k2 | k3 | k4 | y{n+1}")
print("-------------------------------")
#print(Table)
print("")

""" 
-----------ODE To Solve------------------
y' = y / (e^x - 1)
-----------------------------------------
"""
def equation(y,t):
    dydx = y / (math.exp(1)-1)
    return dydx

#plugs numbers into equation
def plugAndChug(x, y):
    yPrime = y / ( math.exp(x) -1 )
    return yPrime


def rungeK(x,y,h):
    
    #figure inital values
    k1 = h * plugAndChug(x,y)
    k2 = h * plugAndChug(x + h / 2, y + k1 / 2)
    k3 = h * plugAndChug(x + h / 2, y + k2 / 2)
    k4 = h * plugAndChug(x + h,y + k3)
    nextY = y + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
    
    
    print("-------------------------------------------------------------------")
    print("------------------Using Runge Kutta Method-------------------------")
    print("Evaluated up to n = 5")
    print("")
    
    #initial values
    n = 1
    nextY = y + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
    print("0", " | ", x, " | ", y, " | ", k1,  " | ", k2, " | ", k3, " | ", k4, " | ", nextY)
    print("")
    #while loop
    while (n < 6):
        #calculate x and y to continue problem
        x = x + h
        y = y + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
        
        #calculate k1,k2,k3 and k4 for next problem
        k1 = h * plugAndChug(x,y)
        k2 = h * plugAndChug(x + h / 2, y + k1 / 2)
        k3 = h * plugAndChug(x + h / 2, y + k2 / 2)
        k4 = h * plugAndChug(x + h,y + k3)
        nextY = y + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
        #print values
        print(n, " | ", x, " | ", y, " | ", k1,  " | ", k2, " | ", k3, " | ", k4, " | ", nextY)
        print("")
        n = n + 1


rungeK(x0,y0,h)


""" Graph ODE using ODEINT """

t = np.linspace(0,2000)
y = odeint(equation,y0,t)

print("-----Data Values Gathered-----")
print("using ODIENT")
print("------------------------------")
print(y)

plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()
