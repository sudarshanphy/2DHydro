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
from mpi4py import MPI

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.family'] = 'DeJavu Serif'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['font.size'] = 18.0

comm = MPI.COMM_WORLD

# Rank of the processor
my_rank = comm.Get_rank()

# No of processors to use
nprocs = comm.Get_size()

# No of files to plot
######################################
nbegin = 0
nend = 2501
nfiles = nend - nbegin + 1

######################################

# Determine the no of files to be provided to each processor
remainder = nfiles % nprocs
quotient = (nfiles - remainder) / nprocs

# Array of number of files for rach processor
files_array = quotient * np.ones(nprocs)

# No of files to be read by each processor
for j in range(remainder):
    files_array[j] = files_array[j] + 1

#def get_blk_info(outdir, basenm, fieldname, num):
#    fname = outdir+basenm+"_0000_%04d.dat"%(num)
#    f = open(fname, "r")
#    lines = f.readlines()[0:7]
#    linem1 = lines[1]
#    line0 = lines[2]
#    line1 = lines[3]
#    line2 = lines[4]
#    
#    time = float(lines[0].split()[2])
#    nx = int(linem1.split()[2])
#    ny = int(line0.split()[2])
#    lnx = int(line1.split()[2])
#    lny = int(line2.split()[2])
#    xblk = int(line1.split()[-1])
#    yblk = int(line2.split()[-1])
#
#    fields = lines[6].split()
#    for i,field in enumerate(fields):
#        if (field == fieldname.lower()):
#            index = i - 1
#
#    print(index, fieldname, time)
#    
#
#    return [nx, ny, lnx, lny, xblk, yblk, index, time]


basenm = "lw_tvd_hllc_x4_y4"
outputdir = "../outputLW/"
fieldname = "dens"
#outputdir = "./"

for nums in range(int(files_array[my_rank])):

    num = my_rank + nbegin #int(sys.argv[1])
    
    fname = outputdir+basenm+"_0000_%04d.dat"%(num)
    f = open(fname, "r")
    lines = f.readlines()[0:7]
    linem1 = lines[1]
    line0 = lines[2]
    line1 = lines[3]
    line2 = lines[4]
    
    time = float(lines[0].split()[2])
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

    print(index, fieldname, time)
    

    #nx, ny, lnx, lny, xblk, yblk, index, time = get_blk_info(outputdir, basenm, fieldname, num)
    
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
        f = open(outputdir+file, "r")
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
    
    
    x_unique = np.unique(x_array)
    y_unique = np.unique(y_array)
    
    X, Y = np.meshgrid(x_unique, y_unique)
    dmin = np.min(data2d)
    dmax = np.max(data2d)
    
    print(dmin, dmax, time)
    dmin = 0.4
    dmax = 1.1
    
    fig = plt.figure(figsize=(12,12)) # , layout='constrained')
    ax = plt.axes()
    plt.title("Liska-Wendroff Implosion Problem (TVD, HLLC) at t= %8.3f"%(time), \
              fontsize=18)
    mesh = plt.pcolormesh(X, Y, np.transpose(data2d), cmap='jet', \
                   vmax=dmax, vmin=dmin) # \
                   #norm=mpl.colors.LogNorm(0.12,1.76))
    
    ax.set_aspect('equal', adjustable='box')
    # this somehow makes the colorbar scale match with the y height
    plt.subplots_adjust(left=0.12, \
                        bottom=0.09, \
                        right = 0.88, \
                        top = 0.94) 
    plt.colorbar(label="Density (g/cc)",fraction=0.046, pad=0.04)
    
    outname = "%s_%s_%04d.png"%(fieldname,basenm,num)
    plt.savefig(outname, dpi=200)   
    print(f"Saved {outname}")
    nbegin = nbegin + nprocs
    plt.clf()
    #plt.show()

MPI.Finalize
