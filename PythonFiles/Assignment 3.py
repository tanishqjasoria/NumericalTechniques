#!/usr/bin/env python
# coding: utf-8

# # Numerical Techniques Laboratory
# ## Assignment 3 | Tanishq Jasoria | 16MA20047

# Solve the following differential equation -
# 
# y''' + 4y" + y' - 6y = 1
# 
# y(0) = y'(0) = 0
# 
# y'(1) = 1

# In[29]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rcParams['figure.figsize'] = [10, 15]


# In[30]:


def BlockTridiagonal(A, B, C, D):
    n = len(D)
    _B = np.zeros(A.shape)
    _C = np.zeros(A.shape)
    _D = np.zeros((n,2,1))
    D_out = np.zeros(_D.shape)
    _C[0] = np.linalg.inv(B[0]).dot(C[0])
    _D[0] = np.linalg.inv(B[0]).dot(D[0])
    
    for i in range(1, n):
        _B[i] = B[i] - A[i].dot(_C[i-1])
        _C[i] = np.linalg.inv(_B[i]) .dot(C[i])
        _D[i] = np.linalg.inv(_B[i]).dot(D[i] - A[i].dot(_D[i-1]))
    D_out[n-1] = np.copy(_D[n-1])
    for i in range(n-2, -1, -1):
        D_out[i] = _D[i] - _C[i].dot(D_out[i+1])
    
    return D_out


# In[31]:


def BVP(x0, xn, h):
    n = int(np.ceil((xn - x0)/h))

    A = np.zeros((n-1, 2, 2))
    B = np.zeros((n-1, 2, 2))
    C = np.zeros((n-1, 2, 2))
    D = np.zeros((n-1, 2, 1))
    for i in range(n-1):
        A[i] = np.array([[-1, -h/2], [0, 1/h**2 - 2/h]])
        B[i] = np.array([[1, -h/2], [-6, -2/h**2 + 1]])
        C[i] = np.array([[0, 0], [0, 1/h**2 + 2/h]])
        D[i] = np.array([[0], [1]])
    
    D[n-2] = D[n-2] - np.array([[0], [1/h**2 + 2/h]])
    X = BlockTridiagonal(A, B, C, D)
    y = X[:, 0]
    
    y = np.reshape(y, n-1)
    print(y.shape)
#     print(y)
#     y = np.append(y, ((1/h - 1/h**2 + 4/h - 3*X[n-2, 1]/h**2)/6))
    print(y)
#     print(A)
#     print(B)
#     print(C)
#     print(D)
    return y


# In[32]:


def func(x0, xn, h = 0.1):
    return np.arange(x0, xn, h)
steps = [0.1, 0.05, 0.02]
colors = ['r', 'g', 'b']
x0 = 0
xn = 1
labels = ["step = 0.1", "step = 0.05", "step = 0.02"]
for step in steps: 
    x_range = func(x0, xn, step)
    print(x_range)
    print("Shape of x_range")
    print(x_range.shape)
    y = BVP(x0, xn, step)
    
    y = np.insert(y, 0, 0)
    print(y.shape)
    print(y)
    plt.xlabel('X')
    plt.ylabel('Y')
    i = steps.index(step)
    plt.plot( x_range, y, colors[i])
    plt.savefig("Plot.png")
    plt.gca().legend(labels)


# In[ ]:




