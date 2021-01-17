import cvxpy as cvx
import numpy as np
import math 

def generate_F_W(n_f, n_w, maxValue=10):
    o_f = np.random.randint(low=0,high = maxValue,size=(n_f))
    i_w = np.random.randint(low=0,high = maxValue,size=(n_w))
    #如果不等，则修改至相等
    sum_o_f = np.sum(o_f)
    sum_i_w = np.sum(i_w)
    sub_value = np.sum(o_f)-np.sum(i_w)
    if sub_value >= 0: #工厂产出量多于仓库存储量
        #修改仓库量：
        i_w[n_w -1] += sub_value 
    else: #工厂产出量较少
        o_f[n_f-1] -= sub_value #负数为减

    return o_f, i_w

#假设仓库点，和存储量
#假设工厂位置点， 和生产量
np.random.seed(1) #固定种子点

#假设10个工厂，前面两个为坐标轴， 最后的量为存储量
n_f = 10
n_w = 11

f = np.random.randint(low =0, high = 20, size=(n_f, 2)) #随机产生坐标位置
w = np.random.randint(low =0, high = 20, size=(n_w, 2))
print('loc of factory:', f)
print('loc of warehouse:', w)
o_f, i_w = generate_F_W(n_f, n_w, maxValue=10)  #随机生成货物量与存储量，单体最大数量
print('output of factories and input of warehouse:',o_f, i_w)

#------------------
# f = np.array([[0,1],[0,0]])
# o_f = np.array([10,10])

# f = np.array([[0,1],[0,0],[3,4]])  #所在位置
# o_f = np.array([10,5,5])

# w = np.array([[3,5],[2,4]])
# i_w = np.array([10,10])
#-----------------------

dis = np.zeros((n_f, n_w))
for i in range(n_f):

    for j in range(n_w):
        dx = math.fabs(f[i][0] - w[j][0])
        dy = math.fabs(f[i][1] - w[j][1])
        dis[i][j] = math.sqrt(dx**2 + dy**2)

print(dis)

Cij = cvx.Variable((n_f, n_w))

#求出相应的工厂和仓库货物量, 要求均大于0，以此为限制条件
constraints = [cvx.sum(Cij, axis=0)==i_w, cvx.sum(Cij, axis=1)==o_f, Cij>=0]  

obj = cvx.Minimize(cvx.sum(cvx.multiply(Cij, dis)))

prob = cvx.Problem(obj, constraints)

prob.solve(solver='CVXOPT', kktsolver=cvx.ROBUST_KKTSOLVER, verbose=True) #求解，并记录返回值


print("status:", prob.status) #求解状态
print("optimal value", prob.value)
print("optimal var: f->w \n", Cij.value)
# cost = cij*dis

#计算距离矩阵；p-q

#2个约束条件

#--------plot figure to show-----------
import matplotlib.pyplot as plt

def get_x_y(arr):
    #获取f\w坐标
    f_point_x = []
    f_point_y = []
    for i in range(len(arr)):
        per_f = arr[i]  #三个位置点都分别取出
        f_point_x.append(per_f[0])
        f_point_y.append(per_f[1])
    return f_point_x, f_point_y

f_point_x, f_point_y = get_x_y(f)
w_point_x, w_point_y = get_x_y(w)

plt.scatter(f_point_x,f_point_y, s = o_f*10+10, marker='o', color='r') #s可以根据货物量来更改大小
plt.scatter(w_point_x,w_point_y, s = i_w*10+10, marker='*', color='b')

#给每个点添加上最大存储量和目前存货量的标志：
for i in range(len(f)): #工厂目前具有的产出量
    #调整显示好位置
    plt.text(f_point_x[i], f_point_y[i]+0.05, o_f[i], ha='center', va= 'bottom',fontsize=11)
for i in range(len(w)):  #显示仓库目前能输出的量
    #调整显示好位置
    plt.text(w_point_x[i], w_point_y[i]+0.05, i_w[i], ha='center', va= 'bottom',fontsize=11)
#标注说明

plt.legend(['Cargo volume in factory','Loading capacity of warehouse'],loc='upper left')
plt.xlabel("East / km")
plt.ylabel("North / km")
plt.show()

###绘画运输路线，加上箭头表示
#propcess the data
transResult = np.array(Cij.value)
print(transResult)
transResult = np.around(transResult, decimals=2)  #指定小数位数
print(transResult, np.shape(transResult))
#获取到相应的工厂输出和仓库获取到的货物量
w_get = np.sum(transResult, axis=0)
f_out = np.sum(transResult, axis=1)
print(w_get, f_out)  
num_factory, num_warehouse = np.shape(transResult)

#剩余量计算，作为第二次画图显示
res_o_f = o_f
res_i_w = i_w
#-------------

for nf in range(num_factory):
    for nw in range(num_warehouse):
        #可以根据坐标来获取哪个工厂运输到哪个仓库
        #xytext:可以看为起始点， xy可以看为终点, 第一个内容为传送的货物量
        if transResult[nf][nw] != 0: #表明有货物运输，则将货物运输情况画出
            plt.annotate(' ', xy=(w_point_x[nw],w_point_y[nw]),xytext=(f_point_x[nf],f_point_y[nf]), ha="center", va = 'center',\
            arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))
            #因为注释文本在箭头起始端，产生覆盖，这里更改文本位置
            #默认标注文本在正中间
            plt.text((w_point_x[nw] + f_point_x[nf])/2, (w_point_y[nw] + f_point_y[nf])/2+0.1, f'{transResult[nf][nw]}') 

            #计算运输量,出库后工厂的存储量下降，入库后仓库的可存储量下降
            res_o_f[nf] -= transResult[nf][nw]
            res_i_w[nw] -= transResult[nf][nw]

#--------------
#此处显示运输完成后的货物量
plt.scatter(f_point_x,f_point_y, s = 20, marker='o', color='r') #s可以根据货物量来更改大小，此处固定大小
plt.scatter(w_point_x,w_point_y, s = 20, marker='*', color='b')

#给每个点添加上最大存储量和目前存货量的标志：
for i in range(len(f)): #工厂目前具有的产出量
    #调整显示好位置
    plt.text(f_point_x[i], f_point_y[i]+0.05, res_o_f[i], ha='center', va= 'bottom',fontsize=11)
for i in range(len(w)):  #显示仓库目前能输出的量
    #调整显示好位置
    plt.text(w_point_x[i], w_point_y[i]+0.05, res_i_w[i], ha='center', va= 'bottom',fontsize=11)
#标注说明

plt.legend(['Cargo volume in factory','Loading capacity of warehouse'],loc='upper left')
plt.xlabel("East / km")
plt.ylabel("North / km")
#----------------

plt.show()

#画出最终结果图
