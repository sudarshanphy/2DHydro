import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.family'] = 'DeJavu Serif'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['font.size'] = 20.0

<<<<<<< HEAD
num = int(sys.argv[1])
#num = 80
fname = "./riemann_test_weno5_%04d.dat"%(num)
fieldname = "dens"
=======
fname = sys.argv[1] #"./kh_test_0010.dat"
fieldname = sys.argv[2] #"dens"
>>>>>>> 39a722a (some custom changes for data)

f = open(fname, "r")

lines = f.readlines()

f.close()

time = float(lines[0].split()[2])
step = float(lines[0].split()[-1])
nx = int(lines[1].split()[2])
dx = float(lines[1].split()[-1])
ny = int(lines[2].split()[2])
dy = float(lines[2].split()[-1])
gamma = float(lines[3].split()[-1])


fields = lines[4].split()
for i,field in enumerate(fields):
    if (field == fieldname.lower()):
        index = i - 1

print(index, fieldname)

#data = np.zeros((nx, ny))
x_array = []
y_array = []
data = []

for line in lines[5:]:
    lst = line.split()
    if (lst[0] != "####"):
        x_array.append(float(lst[0]))
        y_array.append(float(lst[1]))
        data.append(float(lst[index]))
    
    
x_array = np.array(x_array)
y_array = np.array(y_array)
x_array = np.unique(x_array)
y_array = np.unique(y_array)
data = np.array(data).reshape((nx, ny))


X, Y = np.meshgrid(x_array, y_array)
dmin = np.min(data)
dmax = np.max(data)
#print(dmin, dmax)

fig = plt.figure(figsize=(12,12)) # , layout='constrained')
ax = plt.axes()
plt.title("2D Riemann problem (WENO5) at t= %8.3e"%(time))
plt.pcolormesh(X, Y, np.transpose(data), cmap='turbo', \
               #vmax=dmax, vmin=dmin) # \
               norm=mpl.colors.LogNorm(0.12,1.76))


<<<<<<< HEAD
ax.set_aspect('equal', adjustable='box')
# this somehow makes the colorbar scale match with the y height
plt.colorbar(label=fieldname,fraction=0.046, pad=0.04)
=======
plt.figure(figsize=(9, 15))
plt.title("RT problem at t= %8.3e"%(time))
plt.pcolormesh(X, Y, np.transpose(data), cmap='viridis', \
               vmax=dmax, vmin=dmin) # \
               #norm=mpl.colors.LogNorm(dmin,dmax))
plt.colorbar(label=fieldname)
>>>>>>> 39a722a (some custom changes for data)
plt.show()
#plt.savefig("dens_riemann_test_weno3_%04d.png"%(num))   