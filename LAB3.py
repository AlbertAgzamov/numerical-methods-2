import matplotlib.pyplot as plt
import numpy as np
from math import e
from collections import Counter



# Аналитически найденное решение краевой задачи:
def y(x):
    return e ** x + e ** (- x) + 3.1 * x * (x - 1) - 2
    
#q(x) в методе прогонки
def q(x):
    return 3.1 * x * (1 - x) + 8.2

# Правая часть уравнения в краевой задаче:
def f(x, y):
    return y + q(x)

# функция, решающая СЛАУ методом прогонки
def solver(a, b, c, d):
        
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
                    
    xc = bc
    xc[-1] = dc[-1] / bc[-1]

    for il in range(nf - 2, -1, -1):
        xc[il] = (dc[il]-cc[il] * xc[il+ 1]) / bc[il]

    return xc

# Вычисление приближенного решения краевой задачи:
def calculate():

    eps = 0.0001
    yb = e + (1 / e) - 2
    rr = np.arange(0, 1, 0.0001)
    y0 = y(0)
    
    while True:
    
        try:
            N = int(input("Please input natural number N: "))
            if N <= 1: raise ValueError
        except ValueError:
            # Сообщение об ошибке, если N не число или не натуральное число, или веществ. число:
            print("Oops!  That was no valid number.  Try again...")
        else: 
        
            n = N - 1
            
            # h - шаг методов
            h = 1 / n
            h_sq = h ** 2 
            b = -(2 + h_sq)
            
            # Метод прогонки
            if N == 2:
                X_T = [0, 1]
                Y_T = [y0, e + 1 / e - 2]
            elif N == 3:
                X_T = [0, 0.5, 1]
                Y_T = [y0, (q(0.5) * h_sq - 0.7 - e - 1 / e + 2) / b, e + 1 / e - 2]
            else:
                A = [1 for i in range(n - 2)]
                B = [b for i in range(n - 1)]
                C = [1 for i in range(n - 2)]
                x_T = 0
                X_T = [x_T]
                for i in range(n):
                    x_T += h
                    X_T.append(x_T)
                D = []
                i = 1
                for x_T in X_T[1:n]:
                    d = q(x_T) * h_sq
                    if i == n - 1: 
                        d -= e + 1 /e - 2
                    D.append(d)
                    i += 1
                Y_T = [y0]
                Y_T += list(solver(A,B,C,D))
                Y_T += [e + 1 / e - 2]
            
            
            
            # Инициализируем параметр для метода стрельбы (тангенс наклона в начальной точке)
            a = 1
            

            while True:
                
                # Инициализируем начальное значение x (x0)
                x = 0
                
                # Инициализируем значения (y,z)(x0)
                y_S = y0
                z_S = a 
                
                # Инициализируем массивы, куда будем аккумулировать значения y_k
                X = []
                Y_S = []
                X.append(0)
                Y_S.append(y_S)
                
                # Делаем пересчет значений y_S и z_S по методу Рунге-Кутта 4 порядка
                for i in range(n):                  
                    k1 = h * z_S
                    r1 = h * f(x, y_S)
                    
                    k2 = h * (z_S + (1 / 2) * r1)
                    r2 = h * f(x + h / 2, y_S + (1 / 2) * k1)
                    
                    k3 = h * (z_S + (1 / 2) * r2)
                    r3 = h * f(x + h / 2, y_S + (1 / 2) * k2)
                    
                    k4 = h * (z_S + r3)
                    r4 = h * f(x, y_S + k3)
                    
                    x += h
                    
                    
                    y_S += (1 / 6) * k1 + \
                           (1 / 3) * k2 + \
                           (1 / 3) * k3 + \
                           (1 / 6) * k4
                           
                    z_S += (1 / 6) * r1 + \
                           (1 / 3) * r2 + \
                           (1 / 3) * r3 + \
                           (1 / 6) * r4
                    
                    # Добавляем следующие значения в соответствующие массивы значений xk и yk
                    X.append(x)
                    Y_S.append(y_S)
                        
                # Делаем пересчет значений x и y и z и методу трехд. прогонки
                
                if abs(y_S - yb) < eps:
                    # Cюда попадаем, если достигли нужной точности. Нашли решение краевой задачи. Пора строить график.
                    
                    # Cоздаем графики для метода стрельбы, метода прогонки и для точного решения
                    plt.plot(X, Y_S, color='#F01F1F', linestyle='--', marker='.', label='Метод стрельбы')
                    plt.plot(X_T, Y_T, color='#3AF52D', linestyle='--', marker='.', label='Метод прогонки')
                    plt.plot(rr, y(rr), label='Точное решение')
                    
                    # Создаем название для графика
                    plt.title("Сравнение численных методов решения краевой задачи")
                
                    # Даем имена осям
                    plt.xlabel("ось x")
                    plt.ylabel("ось y")
                    
                    # Используем этот метод для удаления лишнего белого/пустого пространства
                    plt.tight_layout()
                    
                    # Стиль графика:
                    plt.style.use('fast')
                    
                    # Добавляем сетку на график, чтоб было удобнее его анализировать
                    plt.grid()
                    
                    # Добавляем легенду
                    plt.legend()
                    
                    # Показываем график
                    plt.show()
                    
                    R = []
                    for i in range(len(Y_T)):       # 1 метод - метод стрельбы, 2 метод - м. прогонки
                        s = abs(Y_S[i] - y(X[i])) # величина отклонения решения 1 метода от точного реш.
                        t = abs(Y_T[i] - y(X[i])) # величина отклонения решения 1 метода от точного реш.
                        if s > t:
                            R.append("T")
                        elif s == t:
                            R.append('E')
                        else: R.append('S')
                    print('In non-boundary values of x:')  # Ниже представлен подсчет количества внутренних узлов, где "побеждает" 1 или второй метод
                    print('Thomas method wins', Counter(R[1:-1])['T'], 'times.') 
                    print('Shooting method wins', Counter(R[1:-1])['S'], 'times.', end='\n\n')
                    break
                
                else:
                    # Cчитаем F(a)
                    Fa = y_S - yb
                    
                    # Считаем F'(a). В коде --> u. 
                    #   Для этого применяем метод Рунге-Кутта 4 порядка
                    u = 0
                    v = 1
                    
                    for j in range(n):
                        
                        q1 = h * v
                        p1 = h * u
                        
                        q2 = h * (v + (1 / 2) * p1)
                        p2 = h * (u + (1 / 2) * q1)
                        
                        q3 = h * (v + (1 / 2) * p2)
                        p3 = h * (u + (1 / 2) * q2)
                        
                        q4 = h * (v + p3)
                        p4 = h * (u + q3) 
                        

                        u += (1 / 6) * q1 + \
                             (1 / 3) * q2 + \
                             (1 / 3) * q3 + \
                             (1 / 6) * q4
                           
                        v += (1 / 6) * p1 + \
                             (1 / 3) * p2 + \
                             (1 / 3) * p3 + \
                             (1 / 6) * p4
                    
                    # Пересчитываем значение параметра a по методу Ньютона-Рафсона
                    a -= Fa / u

calculate()    
            
                