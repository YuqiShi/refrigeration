from CoolProp.CoolProp import PropsSI as CP
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd
import numpy as np
import math

librname = lambda x: 'INCOMP::LiBr[{}]'.format(x)

def h_LiBr(t,x):
    return CP('H','T',t+273.15,'Q',0,librname(x))/1000                      #kJ/kg

def p_LiBr(t,x):
    return CP('P','T',t+273.15,'Q',0,librname(x))/1000                      # kPa
#________________________________________已知x和p计算饱和溶液温度_________________________________________#
def g(t,p,x):
    return p_LiBr(t,x) - p
def t_LiBr(p, x):  
    aa = tc
    bb = tg
    if g(aa,p,x) * g(bb,p,x) >= 0:
        print("NONE")
        return None   
    while (bb - aa) / 2 > 1e-7:
        cc = (aa + bb) / 2
        if g(cc,p,x) == 0:
            return cc
        if g(cc,p,x) * g(aa,p,x) < 0:
            bb = cc
        else:
            aa = cc
    return (aa + bb) / 2
#______________________________________________________________________________________________________#

def s_LiBr(t,x):
    return CP('S','T',t+273.15,'Q',0,librname(x))/1000                      # kJ/kg-K

def f(t,p,x):
    return p_LiBr(t,x) - p
def x_LiBr(t, p):                                                           # kg/kg
    a = 0.2
    b = 0.7
    if f(t,p,a) * f(t,p,b) >= 0:
        return None
    while (b - a) / 2 > 1e-7:
        c = (a + b) / 2
        if f(t,p,c) == 0:
            return c
        if f(t,p,c) * f(t,p,a) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2

def calculate_exit_temperatures(eta, th_in, tc_in, mh, mc, ch, cc):
    # 计算热源和冷源的热容量
    c = min(ch, cc)
    
    # 计算温差
    delta_T = th_in - tc_in
    
    # 计算理论最大传热量
    Q_max = c * delta_T
    
    # 计算实际传热量
    Q = Q_max * eta / 100
    
    # 计算热源出口温度
    th_out = th_in - Q / (mh * ch)
    
    # 计算冷源出口温度
    tc_out = tc_in +(mh * ch / mc / cc) * (th_in - th_out)

    return th_out, tc_out


# External Parameter for a single effect absorption refrigeraiton
tc = 35             # °C
te = 5              # °C
tg = 85             # °C
eta = 0.8           # Heat exchanger efficiency

ph = CP('P','T',tc + 273.15,'Q',0.0,"H2O")/1000
pl = CP('P','T',te + 273.15,'Q',0.0,"H2O")/1000
xh = x_LiBr(tg,ph)
xl = x_LiBr(tc,pl)
m_water = 1
m_weak = xl/(xh-xl)*m_water
m_strong = xh/(xh-xl)*m_water
f = m_weak/m_water

print(xh,xl)
print("Weak Solution: ",m_weak,"Strong Solution: ",m_strong,"Ratio: ",f)


# Initialization
x = np.zeros(10)
t = np.zeros(10)
h = np.zeros(10)
p = np.zeros(10)
s = np.zeros(10)


t[0] = t[7] = tc
t[3] = tg
t[8] = t[9] = te
t[6] = t_LiBr(ph,xl)


p[0] = p[5] = p[8] = p[9] = pl
p[1] = p[2] = p[3] = p[4] = p[6] = p[7] = ph

x[0] = x[1] = x[2] = xl
x[3] = x[4] = x[5] = xh

h[0] = h_LiBr(t[0],xl)

# print(t_LiBr(ph,x[2]))

# Solution Pump
s[0] = s_LiBr(t[1],x[1])
s[1] = s[0]
density1 = CP('D','P',ph*1000,'S',s[1]*1000,librname(xl))
density2 = CP('D','P',pl*1000,'S',s[1]*1000,librname(xl))
w = ph*1000/density1-1000*pl/density2
h[1] = h[0] + w
# print(w, "W/kg",density1,"kg/m^3",density2,"kg/m^3",ph,"kPa")

# Condenser
h[6] = CP('H','T',t[6]+273.15,'P',ph,'H2O')/1000
h[7] = CP('H','T',t[7]+273.15,'Q',0.0,'H2O')/1000
Qc = (h[6] - h[7])*m_water/1000
print(Qc,"kJ:Qc")

# Refrigerant Valve
h[8] = h[7]
Quality = CP('Q','P|twophase',p[8]*1000,'H',h[8]*1000,'Water')
print("Quality:",Quality*100,"%")

# Evaporator
h[9] = CP('H','T',t[9]+273.15,'Q',1.0,'Water')/1000
Qe = (h[9] - h[8])*m_water/1000
print(Qe,"kJ:Qe")

# Solution Heat Exchanger
t[1] = CP('T','P',ph*1000,'H',h[1]*1000,librname(xl))-273.15
CP_ss = CP('C','P',ph*1000,'H',h[1]*1000,librname(xl))/1000
CP_ws = CP('C','P',ph*1000,'H',h[3]*1000,librname(xl))/1000
t[4], t[2] = calculate_exit_temperatures(eta*100, t[3], t[1], m_weak, m_strong, CP_ws, CP_ss)
Q_SHX = (t[2] - t[1]) * CP_ss * m_strong
h[2] = (h[1] * m_strong + Q_SHX)/m_strong
#deltaT = min((1/CP_ss/m_strong),(1/CP_ws/m_weak))
# print(Q_SHX,"kJ")
# print(h[2],"kJ/kg")
# h[2] = CP('H','P',ph*1000,'T',t[2]+273.15,librname(xl))/1000
# print(h[2],"kJ/kg")

# Generator
h[3] = CP('H','T',tg+273.15,'Q',0.0,librname(xh))/1000
Qg = (m_weak * h[3] + m_water * h[6] - m_strong * h[2])/1000
print(Qg,"kJ:Qg")

# Solution Valve
h[4] = CP('H','P',ph*1000,'T',t[4]+273.15,librname(xh))/1000
h[5] = h[4]
#Quality5 = CP('Q','P',p[5]*1000,'H',h[5]*1000,librname(xh))
t[5] = t_LiBr(p[5],x[5])

# Absorber
Qa = (m_water * h[9] + m_weak * h[5] - m_strong * h[0])/1000
print(Qa,"kJ:Qa")


#Entropy
s[2] = s_LiBr(t[2],x[2])
s[3] = s_LiBr(t[3],x[3])
s[4] = s_LiBr(t[4],x[4])
s[5] = s_LiBr(t[5],x[5])

s[6] = CP('S','P',p[5]*1000,'H',h[6]*1000,'Water')/1000
s[7] = CP('S','P|twophase',p[7]*1000,'H',h[7]*1000,'Water')/1000
s[8] = CP('S','P|twophase',p[8]*1000,'H',h[8]*1000,'Water')/1000
s[9] = CP('S','P|twophase',p[9]*1000,'H',h[9]*1000,'Water')/1000

##################  Performance ##################
COP = Qe/(Qg)
COP_e = Qe/(Qg+w*m_strong/1000)
COP_c = (Qc+Qa)/(Qg+w*m_strong/1000)
COP_mechnical = Qe/(w*m_strong/1000)
print("##################################################")
print("###              理论制冷性能：{:.4f}".format(COP),"         ###")
print("###              制冷性能：{:.4f}".format(COP_e),"             ###")
print("###              制热性能：{:.4f}".format(COP_c),"             ###")
print("###              机械制冷性能：{:.4f}".format(COP_mechnical),"        ###")
print("##################################################")
##################################################


data = [t, p, x, h, s]
table = tabulate(data, headers=['0','1', '2', '3','4','5','6','7','8','9'])
print(table)


#########################    State Points Table Ploting     #######################
# 创建一个DataFrame对象，用于存储表格数据
df = pd.DataFrame({
    'StatePoints': ['0','1', '2', '3','4','5','6','7','8','9'],
    'Temperature/°C': t,
    'Pressure/kPa': p,
    'Concentration/100%':x,
    'Enthalpy/kJ·kg-1': h,
    'Entropy/ kJ·kg-1·K-1':s,
})

# 绘制表格
table = plt.table(cellText=df.values, colLabels=df.columns, loc='center')

# 调整表格的位置和大小
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8)
table.scale(3, 3)

# 添加左侧说明文字
#plt.text(0.1, 0.5, '这是左侧说明文字', fontsize=14, ha='center', va='center', rotation='vertical')

# 隐藏坐标轴
plt.axis('off')
# 显示图表
plt.show()