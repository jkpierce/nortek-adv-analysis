import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import adv_analysis as adv
 
fname_z = 'example-data/stage_height.csv'
fname_v = 'example-data/test20220408191222.dat'

dat = adv.collector(fname_z, fname_v,corr_thresh=50,lam_a = 2,lam_v = 1)
stat = adv.binner(dat,fname_v,23)
adv.profile_fitter(dat,stat,fname_v,zmin=0.6,zmax=4.5)
