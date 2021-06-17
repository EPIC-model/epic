from tools.h5_reader import H5Reader
from tools.plot_beautify import *
import matplotlib.colors as cls
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob

filedir='../time_step_tests'
files_ls_rk4=sorted(glob.glob(filedir+'/*ls_rk4*.hdf5'))
files_classic_rk4=sorted(glob.glob(filedir+'/*classic_rk4*.hdf5'))
ls_rk4_ref=filedir+'/tg2d_ls_rk4_0_00078125.hdf5'
classic_rk4_ref=filedir+'/tg2d_classic_rk4_0_00078125.hdf5'

step=4 # Use 0.4 time units, i.e. before first splitting
time_base=0.1
eps=0.01

def plot_error(ref_file,file_list,zorder=0):
    dts=[]
    errors=[]
    h5reader = H5Reader()
    h5reader.open(ref_file)
    ref_positions=h5reader.get_dataset(step, 'position')
    h5reader.close()
    for this_file in file_list:
        h5reader = H5Reader()
        h5reader.open(this_file)
        #num_parcels = h5reader.get_num_parcels(step)
        dt=h5reader.get_step_attribute(step, 'dt')
        #t=h5reader.get_step_attribute(step, 't')
        positions=h5reader.get_dataset(step, 'position')
        #print(this_file,num_parcels,dt,t)
        h5reader.close()
        dts.append(dt[0])
        errors.append(np.sqrt(sum(sum((positions-ref_positions)*(positions-ref_positions)))))
    plt.loglog(dts,errors,'o',zorder=zorder)

plot_error(ls_rk4_ref,files_ls_rk4,zorder=5)
plot_error(classic_rk4_ref,files_classic_rk4,zorder=4)
plot_error(classic_rk4_ref,files_ls_rk4)
plot_error(ls_rk4_ref,files_classic_rk4)
plt.title('error after '+str(time_base*step)+' time units')
plt.loglog([0.1,0.001],[eps*0.1**4,eps*0.001**4])
plt.legend(('ls','classic','ls (classic ref)','classic (ls ref)','4th order slope'))
plt.xlabel('time step length')
plt.ylabel('L2 error wrt reference')
plt.savefig('time-step-checks')
