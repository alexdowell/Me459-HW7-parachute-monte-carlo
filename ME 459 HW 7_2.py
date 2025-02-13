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

 # Drogue to main transition altitude
start_alt = 5000.0 # Starting altitude (Z) location
start_x = 0.0 # Starting X location 
start_y = 0.0 # Starting Y location
#runin_angle = 45*DegtoRad # Heading of aircraft when payload is released
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
        wind_interp_x[i] = np.interp(z[i],wind_data[:,0],wind_data[:,1]) + wind_pdf.rvs() #+(np.random.chisquare(2,len(wind_data[:,1]))) #
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
def box_check(x1,y1,x2,y2,x,y):
    if x > x1: 
        if x < x2: 
            if y > y1: 
                if y < y2:
                    return False #if false its inside the box
                else:
                    return True
            else:
                return True
        else:
            return True
    else:
        return True
#def box_check(x1,y1,x2,y2,x,y):
#    if x>x1 and x < x2 and y>y1 and y<y2:
#        return False #if false its inside the box
#    else:
#        return True
####################################
##### Main Simulation/Program ######
####################################
x_target = 2005
y_target = 1088
number_of_simulations = 5
#z_transition_list = [5000, 4000, 3000, 2000, 1750]
x1 = 2030
x2 = 2130
y1 = 1110
y2 = 1200
keepout_x =[x1,x2]
keepout_y =[y1,y2]
cnt_36 = 0
x_run36 = []
y_run36 = []
distance_run36 = []
cnt_72 = 0
x_run72 = []
y_run72 = []
distance_run72 = []
cnt_108 = 0
x_run108 = []
y_run108 = []
distance_run108 = []
cnt_144 = 0
x_run144 = []
y_run144 = []
distance_run144 = []
cnt_180 = 0
x_run180 = []
y_run180 = []
distance_run180 = []
cnt_216 = 0
x_run216 = []
y_run216 = []
distance_run216 = []
cnt_252 = 0
x_run252 = []
y_run252 = []
distance_run252 = []
cnt_288 = 0
x_run288 = []
y_run288 = []
distance_run288 = []
cnt_324 = 0
x_run324 = []
y_run324 = []
distance_run324 = []
cnt_360 = 0
x_run360 = []
y_run360 = []
distance_run360 = []
ang_i = 360/10
ang_list = []
for ang in range(0,10):
    ang_list.append(ang_i)
    ang_i = ang_i + 360/10
for runin_angle in ang_list:
    for sim in range(number_of_simulations):
        z_transition = 3000
        z_transition_deviation = 75
        z_transition = np.random.normal(z_transition, z_transition_deviation)
        runin_angle1 = runin_angle*DegtoRad
#        runin_angle_dev = 5*DegtoRad
#        runin_angle = np.random.normal(runin_angle,runin_angle_dev)
        runin_speed = np.random.normal(runin_speed,6)
        Cd = np.random.normal(Cd,.1)
        Vx_init = runin_speed*math.cos(runin_angle1) # Take runin angle and speed and convert to x and y initial speeds
        Vy_init = runin_speed*math.sin(runin_angle1)
        
        z,zd, zdd, t = z_sim(start_alt, end_alt, dt ,Sdrogue, Smain, Cd, Mt, z_transition) # Run z (altitude) simulation
        x, xd, xdd, y, yd, ydd = x_y_sim(z, t, dt, start_x, start_y, Vx_init, Vy_init, Cd, Sdrogue, Smain, Mt, z_transition, wind_data) # Run x (downrange) simulation
        if runin_angle == 36.0:
            x_run36.append(x[-1])
            y_run36.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run36.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_36 = cnt_36 +1
        if runin_angle == 72.0:
            x_run72.append(x[-1])
            y_run72.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run72.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_72 = cnt_72 +1
        if runin_angle == 108.0:
            x_run108.append(x[-1])
            y_run108.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run108.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_108 = cnt_108 +1
        if runin_angle == 144.0:
            x_run144.append(x[-1])
            y_run144.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run144.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_144 = cnt_144 +1
        if runin_angle == 180.0:
            x_run180.append(x[-1])
            y_run180.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run180.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_180 = cnt_180 +1
        if runin_angle == 216.0:
            x_run216.append(x[-1])
            y_run216.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run216.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_216 = cnt_216 +1
        if runin_angle == 252.0:
            x_run252.append(x[-1])
            y_run252.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run252.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_252 = cnt_252 +1
        if runin_angle == 288.0:
            x_run288.append(x[-1])
            y_run288.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run288.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_288 = cnt_288 +1
        if runin_angle == 324.0:
            x_run324.append(x[-1])
            y_run324.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run324.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_324 = cnt_324 +1
        if runin_angle == 360.0:
            x_run360.append(x[-1])
            y_run360.append(y[-1])
            distance= np.sqrt((x[-1]-x_target)*(x[-1]-x_target)+(y[-1]-y_target)*(y[-1]-y_target))
            distance_run360.append(distance)
            if box_check(x1,y1,x2,y2,x[-1],y[-1]) == False:
                cnt_360 = cnt_360 +1
mean_distance_from_target_36 = sum(distance_run36)/len(distance_run36)
mean_distance_from_target_72 = sum(distance_run72)/len(distance_run72)
mean_distance_from_target_108 = sum(distance_run108)/len(distance_run108)
mean_distance_from_target_144 = sum(distance_run144)/len(distance_run144)
mean_distance_from_target_180 = sum(distance_run180)/len(distance_run180)
mean_distance_from_target_216 = sum(distance_run216)/len(distance_run216)
mean_distance_from_target_252 = sum(distance_run252)/len(distance_run252)
mean_distance_from_target_288 = sum(distance_run288)/len(distance_run288)
mean_distance_from_target_324 = sum(distance_run324)/len(distance_run324)
mean_distance_from_target_360 = sum(distance_run360)/len(distance_run360)
standard_deviation_distance_from_target_36 = statistics.pstdev(distance_run36)
standard_deviation_distance_from_target_72 = statistics.pstdev(distance_run72)
standard_deviation_distance_from_target_108 = statistics.pstdev(distance_run108)
standard_deviation_distance_from_target_144 = statistics.pstdev(distance_run144)
standard_deviation_distance_from_target_180 = statistics.pstdev(distance_run180) 
standard_deviation_distance_from_target_216 = statistics.pstdev(distance_run216)
standard_deviation_distance_from_target_252 = statistics.pstdev(distance_run252)
standard_deviation_distance_from_target_288 = statistics.pstdev(distance_run288)
standard_deviation_distance_from_target_324 = statistics.pstdev(distance_run324)
standard_deviation_distance_from_target_360 = statistics.pstdev(distance_run360)                                             
percent_outside_box_36 = (len(distance_run36)-cnt_36) / len(distance_run36)
percent_outside_box_72 = (len(distance_run72)-cnt_72) / len(distance_run72)
percent_outside_box_108 = (len(distance_run108)-cnt_108) / len(distance_run108)
percent_outside_box_144 = (len(distance_run144)-cnt_144) / len(distance_run144)
percent_outside_box_180 = (len(distance_run180)-cnt_180) / len(distance_run180)
percent_outside_box_216 = (len(distance_run216)-cnt_216) / len(distance_run216)
percent_outside_box_252 = (len(distance_run252)-cnt_252) / len(distance_run252)
percent_outside_box_288 = (len(distance_run288)-cnt_288) / len(distance_run288)
percent_outside_box_324 = (len(distance_run324)-cnt_324) / len(distance_run324)
percent_outside_box_360 = (len(distance_run360)-cnt_360) / len(distance_run360)
table = [['Run in Angle(deg)', 'Mean Distance from Target(m)', 'Standard Deviation from Target(m)', '% Outside No-Go Zone'], ['36', mean_distance_from_target_36, standard_deviation_distance_from_target_36,percent_outside_box_36],['72', mean_distance_from_target_72, standard_deviation_distance_from_target_72,percent_outside_box_72],['108', mean_distance_from_target_108, standard_deviation_distance_from_target_108,percent_outside_box_108],['144', mean_distance_from_target_144, standard_deviation_distance_from_target_144,percent_outside_box_144],['180', mean_distance_from_target_180, standard_deviation_distance_from_target_180,percent_outside_box_180],['216', mean_distance_from_target_216, standard_deviation_distance_from_target_216,percent_outside_box_216],['252', mean_distance_from_target_252, standard_deviation_distance_from_target_252,percent_outside_box_252],['288', mean_distance_from_target_288, standard_deviation_distance_from_target_288,percent_outside_box_288],['324', mean_distance_from_target_324, standard_deviation_distance_from_target_324,percent_outside_box_324],['360', mean_distance_from_target_360, standard_deviation_distance_from_target_360,percent_outside_box_360]]
print(tabulate(table, headers='firstrow'))
plt_36 = plt.scatter(x_run36, y_run36, s = 1, c ="blue")
plt_72 = plt.scatter(x_run72, y_run72, s = 1, c ="pink")
plt_108 = plt.scatter(x_run108, y_run108, s = 1, c ="orange")
plt_144 = plt.scatter(x_run144, y_run144, s = 1, c ="green")
plt_180 = plt.scatter(x_run180, y_run180, s = 1, c ="yellow")
plt_216 = plt.scatter(x_run216, y_run216, s = 1, c ="red")
plt_252 = plt.scatter(x_run252, y_run252, s = 1, c ="black")
plt_288 = plt.scatter(x_run288, y_run288, s = 1, c ="purple")
plt_324 = plt.scatter(x_run324, y_run324, s = 1, c ="brown")
plt_360 = plt.scatter(x_run360, y_run360, s = 1, c ="cyan")
plt.plot([keepout_x[0],keepout_x[1], keepout_x[1], keepout_x[0], keepout_x[0]], [keepout_y[0], keepout_y[0], keepout_y[1], keepout_y[1], keepout_y[0]], 'r-') 
plt.legend((plt_36, plt_72, plt_108, plt_144, plt_180,plt_216,plt_252,plt_288,plt_324,plt_360,),
           ('36', '72', '108', '144', '180', '216', '252', '288', '252', '288', '324', '360'),
           scatterpoints=1,
           loc='upper left',ncol=5,
           fontsize=8)
plt.xlabel("x-distance (m)")
plt.ylabel("y-distance (m)")
plt.title("Landing Locations")
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






