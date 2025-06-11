'''
03/09/2025
Sudarshan Neopane
This file is to read in data from a parallel run.
Each block write's it's own output file.
For a block with rank 4, the outfile file at initilization will be
filename: basenm_0004(block rank)_0000.dat
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.family'] = 'DeJavu Serif'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['font.size'] = 20.0

def get_blk_info(basenm, fieldname, num):
    fname = "../output/"+basenm+"_0000_%04d.dat"%(num)
    f = open(fname, "r")
    lines = f.readlines()[0:7]
    linem1 = lines[1]
    line0 = lines[2]
    line1 = lines[3]
    line2 = lines[4]
    
    nx = int(linem1.split()[2])
    ny = int(line0.split()[2])
    lnx = int(line1.split()[2])
    lny = int(line2.split()[2])
    xblk = int(line1.split()[-1])
    yblk = int(line2.split()[-1])

    fields = lines[6].split()
    for i,field in enumerate(fields):
        if (field == fieldname.lower()):
            index = i - 1

    print(index, fieldname)
    

    return [nx, ny, lnx, lny, xblk, yblk, index]


def get_data(basenm, num, fieldname):

    nx, ny, lnx, lny, xblk, yblk, index = get_blk_info(basenm, fieldname, num)
    # get all the file names
    fnames = []
    count = 0
    for j in range(yblk):
        for i in range(xblk):
            fnames.append(basenm+"_%04d_%04d.dat"%(count,num))
            count = count + 1

    x_array = []
    y_array = []

    data2d = np.zeros((nx,ny))

    for kk,file in enumerate(fnames):
        f = open("../output/"+file, "r")
        #print(f)
        lines = f.readlines()
        
        f.close()
        data = []
        xrank = int(kk%xblk)
        yrank = int(np.floor(kk/xblk))
        ilo = int(xrank * lnx)
        ihi = int((xrank + 1)*lnx - 1)
        jlo = int(yrank * lny)
        jhi = int((yrank + 1)*lny - 1)
        #print(xrank, yrank)
        #print(ilo, ihi, jlo, jhi)
        for line in lines[7:]:
            lst = line.split()
            if (lst[0] != "####"):
                x_array.append(float(lst[0]))
                y_array.append(float(lst[1]))
                try:
                    data.append(float(lst[index]))
                except:
                    data.append(float(0.0))
        data2d[ilo:ihi+1,jlo:jhi+1] = np.array(data).reshape((lnx,lny))

    return [x_array, y_array, data2d]


basenm1 = "sedov_newrecon_test_x2_y2"
basenm2 = "sedov_test_x2_y2"
num = 11
fieldname = "dens"

x_array, y_array, data2d = get_data(basenm1, num, fieldname)
x_array2, y_array2, data2d2 = get_data(basenm2, num, fieldname)

x_unique = np.unique(x_array)
y_unique = np.unique(y_array)

X, Y = np.meshgrid(x_unique, y_unique)

diffdata = (data2d - data2d2)

dmin = np.min(diffdata)
dmax = np.max(diffdata)

print(dmin, dmax)


fig = plt.figure(num = basenm1, figsize=(12,12)) # , layout='constrained')
ax = plt.axes()
#plt.title("2D Riemann problem (WENO5) at t= %8.3e"%(time))
mesh = plt.pcolormesh(X, Y, np.transpose(diffdata), cmap='jet', \
               vmax=dmax, vmin=dmin) # \
               #norm=mpl.colors.LogNorm(0.12,1.76))

ax.set_aspect('equal', adjustable='box')
# this somehow makes the colorbar scale match with the y height
plt.colorbar(label=fieldname,fraction=0.046, pad=0.04)


#plt.savefig("pres_otmhd_x2_2_%04d.png"%(num), dpi=144)   
plt.show()
