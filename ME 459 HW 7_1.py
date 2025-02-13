# Uncontrolled Parachute Simulation v1 #
# For use in Robotics and Unmanned Systems homework #
#  By Travis Fields - March 30, 2020   #

from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import math
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import chi2
import statistics
from tabulate import tabulate

wind_data = genfromtxt('Wind_DataFile.csv', delimiter=',')
wind_pdf = chi2(2,0.1)
# Column 1 = altitude
# Column 2 = wind Vx
# Column 3 = wind Vy

DegtoRad = 3.14/180 # Convert degrees to radians

# Parachute/Payload Parameters #
Sdrogue = 0.3 # Drogue parachute reference area, m^2
Smain = 2.0 # Main parachute reference area, m^2
Cd = 0.7 # Parachute drag coefficient
Mp = 0.5 # Mass of parachute, kg
Ms = 5.0 # Payload mass, kg
Mt = Ms + Mp # Total mass, kg
# ---------------------------- #

# Simulation Parameters #
z_transition_list = [5000, 4000, 3000, 2000, 1750]
 # Drogue to main transition altitude
start_alt = 7000.0 # Starting altitude (Z) location
start_x = 0.0 # Starting X location 
start_y = 0.0 # Starting Y location
runin_angle = 45*DegtoRad # Heading of aircraft when payload is released
runin_speed = 75 # Speed of aircraft [m/s] when payload is released
end_alt = 1500.0 # Ground altitude
dt = 0.05 # Simulation time step, s (if too large simulation goes unstable)
# --------------------- #



def z_sim(z_init, z_end, dt, Sd, Sm, Cd, M, z_trans):
    g = 9.81 # SI units
    z = np.ones(1) # Initializing altitude array, will append as sim runs
    zd = np.ones(1) # Init. z dot array
    zdd = np.ones(1) # Init. z double dot array
    t = np.ones(1) # Initializing time array, will append as sim runs
    t[0] = 0
    z[0] = z_init # Start sim at release altitude
    zd[0] = -0.001 # Starting with 0 descent velocity (may have forward velocity though)
    zdd[0] = 0 # 0 descent acceleration for first time step
	
    i = 0
    while (z[i] > z_end): # Looping until altitude is below the ground
        i = i + 1
        if (z[i-1]>z_trans): # If still above transition altitude
            S = Sd # use drogue parachute size
        else:
            S = Sm # otherwise, use main parachute size
            
        zdd = np.append(zdd, -g - np.sign(zd[i-1])*0.5/M*rho(z[i-1])*Cd*S*zd[i-1]**2)
        # z accel is weight and drag only (no other external force/propulsion)
        
        zd = np.append(zd, zd[i-1] + zdd[i]*dt)
        # z velocity is simple kinematic integration, vel = vel_prev + accel*dt
        
        z = np.append(z, z[i-1] + zd[i]*dt + 0.5*zdd[i]*dt**2)
        # altitude is simple kinematic integration
        
        t = np.append(t, t[i-1] + dt) # Simple, but sticking in here for convenience
    return z, zd, zdd, t


def x_y_sim(z, t, dt, x0, y0, Vx0, Vy0, Cd, Sd, Sm, M, z_trans, wind):
    x = np.ones(1) # Init X position array
    xd = np.ones(1) # Init X velocity array
    xdd = np.ones(1) # Init X acceleration array
    y = np.ones(1) # Init Y position array
    yd = np.ones(1) # Init Y velocity array
    ydd = np.ones(1) # Init Y acceleration array
    
    # Setting initial conditions
    x[0] = x0 # Initial position from aircraft deployment (0,0 is "perfect")
    xd[0] = Vx0 # Initial velocity from aircraft deployment
    xdd[0] = 0 
    y[0] = y0 # Initial position from aircraft deployment (0,0 is "perfect")
    yd[0] = Vy0 # Initial velocity from aircraft deployment
    ydd[0] = 0
    
    wind_interp_x = np.ones(len(z)) # Init wind interpolation (need wind at z(k) altitude, not altitudes given in file)
    wind_interp_x[0] = 0.0
    wind_interp_y = np.ones(len(z))
    wind_interp_y[0] = 0.0
    
    for i in range(1,len(z)):
#        wind_interp_x[i] = np.interp(z[i],wind[:,0], wind[:,1])  #add wind noise here# Does not have any randomness added to interpolation
#        wind_interp_y[i] = np.interp(z[i], wind[:,0], wind[:,2]) 
        # this linear interpolation uses the wind data (frrom ) and the previous altitude (z[i-1]) to estimate what the wind is at z[i-1]
        wind_interp_x[i] = np.interp(z[i],wind_data[:,0],wind_data[:,1]) + wind_pdf.rvs()#+(np.random.chisquare(2,len(wind_data[:,1]))) #
        wind_interp_y[i] = np.interp(z[i], wind_data[:,0],wind_data[:,2]) + wind_pdf.rvs() #+(np.random.chisquare(2,len(wind_data[:,2])))        
        xdd = np.append(xdd,0) # Assuming zero horizontal acceleration
        ydd = np.append(ydd,0)
        
        xd = np.append(xd, wind_interp_x[i]) # Assumes traveling with wind
        yd = np.append(yd, wind_interp_y[i])
        
        x = np.append(x, x[i-1] + xd[i]*dt + 0.5*xdd[i]*dt**2)
        y = np.append(y, y[i-1] + yd[i]*dt + 0.5*ydd[i]*dt**2)
        # position is simply kinematic integration 
    return x, xd, xdd, y, yd, ydd
	
	
def rho(alt):
    C1 = -3.9142e-14
    C2 = 3.6272e-9
    C3 = -1.1357e-4
    C4 = 1.2204
    rho_poly = C1*alt**3 + C2*alt**2 + C3*alt + C4
    return rho_poly

####################################
##### Main Simulation/Program ######
####################################
x_target = 3000
y_target = 2000
number_of_simulations = 100
z_transition_list = [5000, 4000, 3000, 2000, 1750]
x_at_transition_5000 = []
y_at_transition_5000 = []
distance_at_transition_5000 = []
within_tolerance_distance_at_transition_5000 = []
x_at_transition_4000 = []
y_at_transition_4000 = []
distance_at_transition_4000 = []
within_tolerance_distance_at_transition_4000 = []
x_at_transition_3000 = []
y_at_transition_3000 = []
distance_at_transition_3000 = []
within_tolerance_distance_at_transition_3000 = []
x_at_transition_2000 = []
y_at_transition_2000 = []
distance_at_transition_2000 = []
within_tolerance_distance_at_transition_2000 = []
x_at_transition_1750 = []
y_at_transition_1750 = []
distance_at_transition_1750 = []
within_tolerance_distance_at_transition_1750 = []
for z_transition in z_transition_list:
    for sim in range(number_of_simulations):
        
        runin_angle_dev = 5*DegtoRad
        runin_angle = np.random.normal(runin_angle,runin_angle_dev)
        runin_speed = np.random.normal(runin_speed,4)
        Cd = np.random.normal(Cd,.1)
        Vx_init = runin_speed*math.cos(runin_angle) # Take runin angle and speed and convert to x and y initial speeds
        Vy_init = runin_speed*math.sin(runin_angle)
        
        z,zd, zdd, t = z_sim(start_alt, end_alt, dt ,Sdrogue, Smain, Cd, Mt, z_transition) # Run z (altitude) simulation
        x, xd, xdd, y, yd, ydd = x_y_sim(z, t, dt, start_x, start_y, Vx_init, Vy_init, Cd, Sdrogue, Smain, Mt, z_transition, wind_data) # Run x (downrange) simulation
        if z_transition == 5000:
            x_at_transition_5000.append(x[-1])
            y_at_transition_5000.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_at_transition_5000.append(distance)
            if distance <= 100:
                within_tolerance_distance_at_transition_5000.append(distance)
        if z_transition == 4000:
            x_at_transition_4000.append(x[-1])
            y_at_transition_4000.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_at_transition_4000.append(distance)
            if distance <= 100:
                within_tolerance_distance_at_transition_4000.append(distance)
        if z_transition == 3000:
            x_at_transition_3000.append(x[-1])
            y_at_transition_3000.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_at_transition_3000.append(distance)
            if distance <= 100:
                within_tolerance_distance_at_transition_3000.append(distance)
        if z_transition == 2000:
            x_at_transition_2000.append(x[-1])
            y_at_transition_2000.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_at_transition_2000.append(distance)
            if distance <= 100:
                within_tolerance_distance_at_transition_2000.append(distance)
        if z_transition == 1750:
            x_at_transition_1750.append(x[-1])
            y_at_transition_1750.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_at_transition_1750.append(distance)
            if distance <= 100:
                within_tolerance_distance_at_transition_1750.append(distance)
mean_distance_from_target_5000 = sum(distance_at_transition_5000)/len(distance_at_transition_5000)
mean_distance_from_target_4000 = sum(distance_at_transition_4000)/len(distance_at_transition_4000)
mean_distance_from_target_3000 = sum(distance_at_transition_3000)/len(distance_at_transition_3000)
mean_distance_from_target_2000 = sum(distance_at_transition_2000)/len(distance_at_transition_2000)
mean_distance_from_target_1750 = sum(distance_at_transition_1750)/len(distance_at_transition_1750)
standard_deviation_distance_from_target_5000 = statistics.pstdev(distance_at_transition_5000)
standard_deviation_distance_from_target_4000 = statistics.pstdev(distance_at_transition_4000)
standard_deviation_distance_from_target_3000 = statistics.pstdev(distance_at_transition_3000)
standard_deviation_distance_from_target_2000 = statistics.pstdev(distance_at_transition_2000)
standard_deviation_distance_from_target_1750 = statistics.pstdev(distance_at_transition_1750)                                             
percent_within_tolerance_5000 = len(within_tolerance_distance_at_transition_5000) / len(distance_at_transition_5000)
percent_within_tolerance_4000 = len(within_tolerance_distance_at_transition_4000) / len(distance_at_transition_4000)
percent_within_tolerance_3000 = len(within_tolerance_distance_at_transition_3000) / len(distance_at_transition_3000)
percent_within_tolerance_2000 = len(within_tolerance_distance_at_transition_2000) / len(distance_at_transition_2000)
percent_within_tolerance_1750 = len(within_tolerance_distance_at_transition_1750) / len(distance_at_transition_1750)
table = [['Transition Altitude', 'Mean Distance from Target', 'Sandard Deviation Distance', '% in Bounds'], ['5000', mean_distance_from_target_5000, standard_deviation_distance_from_target_5000,percent_within_tolerance_5000],['4000', mean_distance_from_target_4000, standard_deviation_distance_from_target_4000,percent_within_tolerance_4000],['3000', mean_distance_from_target_3000, standard_deviation_distance_from_target_3000,percent_within_tolerance_3000],['2000', mean_distance_from_target_2000, standard_deviation_distance_from_target_2000,percent_within_tolerance_2000],['1750', mean_distance_from_target_1750, standard_deviation_distance_from_target_1750,percent_within_tolerance_1750]]
print(tabulate(table, headers='firstrow'))
plt_5000 = plt.scatter(x_at_transition_5000, y_at_transition_5000, s = 1, c ="blue")
plt_4000 = plt.scatter(x_at_transition_4000, y_at_transition_4000, s = 1, c ="pink")
plt_3000 = plt.scatter(x_at_transition_3000, y_at_transition_3000, s = 1, c ="orange")
plt_2000 = plt.scatter(x_at_transition_2000, y_at_transition_2000, s = 1, c ="green")
plt_1750 = plt.scatter(x_at_transition_1750, y_at_transition_1750, s = 1, c ="yellow")
plt.legend((plt_5000, plt_4000, plt_3000, plt_2000, plt_1750),
           ('5000', '4000', '3000', '2000', '1750'),
           scatterpoints=1,
           loc='upper left',ncol=5,
           fontsize=8)
plt.xlabel("x-distance (m)")
plt.ylabel("y-distance (m)")
plt.title("Landing Location")
plt.show()
#plt.plot(x_at_transition_5000,y_at_transition_5000,'c', label='5000')
#plt.plot(x_at_transition_4000,y_at_transition_4000,'k', label='4000')
#plt.plot(x_at_transition_3000,y_at_transition_3000,'m', label='3000')
#plt.plot(x_at_transition_2000,y_at_transition_2000,'b', label='2000')
#plt.plot(x_at_transition_1750,y_at_transition_1750,'g', label='1750')       
#        fig = plt.figure(1)
#        ax = fig.add_subplot(111,projection='3d')
#        ax.plot3D(x, y, z)
#        ax.set_xlabel('X Distance, m')
#        ax.set_ylabel('Y Distance, m')
#        ax.set_zlabel('Altitude, m')
        
#        fig = plt.figure(2)
#        plt.plot(t,zd)
#        plt.xlabel('Time, s')
#        plt.ylabel('Descent Velocity, m/s')
#        
#        fig = plt.figure(3)
#        plt.plot(t,yd)
#        plt.xlabel('Time, s')
#        plt.ylabel('Y Velocity, m/s')
#        
#        fig = plt.figure(4)
#        plt.plot(t,xd)
#        plt.xlabel('Time, s')
#        plt.ylabel('X Velocity, m/s')




