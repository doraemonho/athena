# Regression test based on the decaying Alfven wave due to viscosity
# and resistivity. The decay rate is fit and then compared with 
# analytic solution

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data
from scipy.optimize import curve_fit

def prepare(**kwargs):
  athena.configure('b',
      prob='linear_wave',
      flux='hlld',
      eos='adiabatic') #isothermal')
  athena.make()

def run(**kwargs):
  arguments0 = ['output1/file_type=hst','output1/dt=0.01',
      'output2/file_type=vtk','output2/variable=prim','output2/dt=0.03',
      #'time/cfl_number=0.3','time/tlim=3.0','time/nlim=-1',
      'time/cfl_number=0.3','time/tlim=6.0','time/nlim=-1',
      'mesh/nx1=64','mesh/x1min=0.0','mesh/x1max=3.0','mesh/ix1_bc=periodic','mesh/ox1_bc=periodic',
      'mesh/nx2=32','mesh/x2min=0.0','mesh/x2max=1.5','mesh/ix2_bc=periodic','mesh/ox2_bc=periodic',
      'mesh/nx3=32','mesh/x3min=0.0','mesh/x3max=1.5','mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
      'mesh/num_threads=1','mesh/refinement=none',
      'meshblock/nx1=64','meshblock/nx2=32','meshblock/nx3=32',
      'hydro/iso_sound_speed=1.0',
      #'problem/compute_error=false','problem/wave_flag=1',
      'problem/compute_error=false','problem/wave_flag=0',
      'problem/amp=1.0e-4','problem/vflow=0.0',
      #'problem/nuiso=0.02','problem/coef_o=0.01']
      'problem/nuiso=0.01','problem/coef_o=0.00','problem/kiso=0.01']
  arguments = arguments0+['job/problem_id=DecayAlfven0']
  athena.run('mhd/athinput.linear_wave3d', arguments)

  # operator split viscosity only
  arguments = arguments0+['problem/operator_split_viscosity=1',
                          'job/problem_id=DecayAlfven1']
  #athena.run('mhd/athinput.linear_wave3d', arguments)

  # operator split Ohmic resistivity only
  arguments = arguments0+['problem/operator_split_field_diffusion=1',
                          'job/problem_id=DecayAlfven2']
  #athena.run('mhd/athinput.linear_wave3d', arguments)

  # operator split for both viscosity and Ohmic resistivity
  arguments = arguments0+['problem/operator_split_field_diffusion=1',
                          'problem/operator_split_viscosity=1',
                          'job/problem_id=DecayAlfven3']
  #athena.run('mhd/athinput.linear_wave3d', arguments)

def analyze():
  ksqr = (2.0*np.pi)**2
  rate = 0.5*(0.02+0.01)*ksqr #(nu+eta)*k^2/2 decay rate of Alfven wave
  rate = 0.5*(0.01*4.0/3.0+0.01*4.0/15.0)*ksqr #(4nu/3+(gamma-1)^2/gamma*kappa)*k^2/2 decay rate of sound wave with thermal conduction

  basename='bin/DecayAlfven0.block0.out2.'
  nframe = 101
  nframe = 201
  dumprate = 0.03
  max_vy = np.zeros(nframe)
  tt     = np.zeros(nframe)
  for i in range(nframe):
    x1f,x2f,x3f,data = athena_read.vtk(basename+str(i).zfill(5)+'.vtk')
    max_vy[i] = np.max(data['vel'][...,1])
    tt[i]     = i*dumprate

  # estimate the decay rate from simulation
  #def func(x,a,b,c):
  #    return a*np.exp(-b*x)+c
  #popt,pcov = curve_fit(func,tt,max_vy)
  #new_rate = popt[1]
  #print '[Decaying Linear Wave-3D]: Ref(decay_rate) = ',rate
  #print '[Decaying Linear Wave-3D]: New(decay_rate) = ',new_rate
  #print 'optimal parameter values: (a,b,c) = ',popt[0],popt[1],popt[2] 

  #def func(x,b):
  #    return max_vy[0]*np.exp(-b*x)
  def func(x,a,b):
      return a*np.exp(-b*x)
  popt,pcov = curve_fit(func,tt,max_vy)
  new_rate = popt[1]
  print '[Decaying Linear Wave-3D]: Ref(decay_rate) = ',rate
  print '[Decaying Linear Wave-3D]: New(decay_rate) = ',new_rate

  flag = True
  error_rel = np.fabs(rate/new_rate-1.0)
  if error_rel > 0.1:
    print('[Decaying Linear Wave-3D]: decay rate is off by 10 percent')
    flag = False
  else:
    print('[Decaying Linear Wave-3D]: decay rate is within 10 percent precision')

  return flag
