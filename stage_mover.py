from zaber_motion import Library, LogOutputMode
from zaber_motion.binary import Connection, CommandCode
from zaber_motion import Units
from time import sleep
import time
import numpy as np
from datetime import datetime
import sys


Library.enable_device_db_store()

zmin = 920
zmax = 968
#dz = 0.5
N = 50
Z = zmax - np.geomspace(1e-1,zmax-zmin,N) # this is geometric with closest spacing near the bed.
#Z = np.arange(zmin,zmax+dz,dz) # set of positions
#print(Z)

O = [] # set of times at which the stage moved in s
fname = sys.argv[1]+'/stage_height.csv'
twait= 10.0 # time in seconds to wait at each position

# give an expectation of how long.
# then require the user to start



def timenow():
    now = datetime.now()
    return (now - now.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()

with Connection.open_serial_port("COM1") as connection: # open the serial port
    device = connection.detect_devices()[0]	# here is the device
    with open(fname,"w") as file: # open the file to write the timeseries of position vs time
        file.write("time_s,position_mm\n") # set up the header to the timeseries file
        for i,z in enumerate(Z):
            z1 = device.move_absolute(z, Units.LENGTH_MILLIMETRES) # set the initial position
            if i==0:
                #print((zmax-zmin)//dz*twait/60, ' minutes required.')
                print(len(Z)*twait/60, 'minutes required. ')
                input("Start the ADV measurement, then press Enter to continue.")
                print('Measurement started.')

            file.write("{},{}\n".format(timenow(),z1))    
            time.sleep(twait)
            
            if i==len(Z)-1:
                input("Stepping finished. End ADV measurement, then press enter to conclude.")
            file.write("{},{}\n".format(timenow(),z1))
            print(timenow(),z1)
            print('step {} of {}'.format(i+1,len(Z)))
            
            
        z1 = device.move_absolute(850, Units.LENGTH_MILLIMETRES)            
            
