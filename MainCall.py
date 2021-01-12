# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 18:54:39 2020

@author: mstew
"""
##                  Final Version 

#





 ##                                 Import Section
import datetime as dt
import numpy as np
import math
from scipy.spatial.transform import Rotation as R
import csv
import os




 ##                   Satellite Position,Velocity Calculator Functions



    # In[]

#Converts time in TLE format to Datetime Format
def refepoch_to_dt(refepoch):
    global dfrac
    
    Epochyrday = dt.datetime.strptime((refepoch[0:5]),'%y%j')
    dfrac = np.modf(np.float(refepoch))[0]
    
    dfracdt = dt.timedelta(microseconds=np.int(dfrac*24*3600*10**6))
    Epochdt = Epochyrday + dfracdt
    
    #Epochdt=Epochdt.replace(microsecond=int(str(Epochdt.microsecond)[0:3]+"0000"))
    #Coverts to STK accuracy
    #This doesn't make a different in Mean ANomaly Calcualtions
    
    
    
    return Epochdt

    # In[]
# Finds the Day of Year from the Year Month and Day parameters
def doy(YR,MO,D):
    if len(str(MO))==1:
        MO='0'+ str(MO)
        #makes sure Month is represented in 2 digits
    if len(str(D))==1:
        D='0'+ str(D)
        #makes sure Day is in two digit format that is readable by strptime
    String=str(YR)+str(MO)+str(D)
    Time = dt.datetime.strptime(String,'%Y%m%d')
    #Converts First to dt
    DOY=dt.datetime.strftime(Time,'%j')
    #Converts from dt to day of the year
    return DOY
    
    # In[]
# This creates a list of datetime objectas in the time_list
def referenceepoch_propagate(TrackingData):
    
    
    
    starttime=dt.datetime.strptime(TrackingData.starttime,'%Y-%m-%d-%H:%M:%S')
    endtime=dt.datetime.strptime(TrackingData.endtime,'%Y-%m-%d-%H:%M:%S')
    timestep=float(TrackingData.timestep)
    Iterations=(endtime-starttime).total_seconds()/float(TrackingData.timestep)
    time=starttime
    global time_list
    time_list=[]
    
    time_list.append(time)
    for i in range(0,int(Iterations-1)):
        time=time+dt.timedelta(seconds=timestep)

        
        time_list.append(time)
        
        
    return time_list

    
    # In[]
def THETAN(time_array):
    
    
    #Input is a refepoch array this dsoesn't make sense as Tracking Data contains an easily parsible 
    
    
    GMST_list=[]
    
    
    J2000=dt.datetime.strptime('2000-01-01 00:00:00','%Y-%m-%d %H:%M:%S')
    global Starttime
    Starttime=dt.datetime.strptime(Tracking.starttime,'%Y-%m-%d-%H:%M:%S')
    
    

    
    

    for i in range(0,len(time_array)):

        Midtime=time_array[i].replace(hour=0,minute=0,second=0)
        D_u=(Midtime-J2000).days
        T_u=D_u/36525
        GMST_00=(99.9677947+36000.7700631*T_u+0.00038793*T_u**2-2.6e-8*T_u**3)%360
        r=1.002737909350795+5.9006e-11*(T_u)-5.9e-15*(T_u)**2
        #Creates T mid for Observation Day
        #Notice how we replace hour,min and sec to 0. This makes the time midnight!
        del_sec=(time_array[i]-Midtime).total_seconds()
        

        

        
        GMST_t=(GMST_00+360*r/86400*(del_sec))%360
        
        GMST_list.append(GMST_t)
    
    
    
        
    return GMST_list    

    # In[]
def mean_anomaly_motion(time,ts_sat_epoch,M0_mean_anomaly,n_mean_motion, \
                        n_dot_mean_motion,n_2dots_mean_motion):
    
    
    #time- datetime object
    # n_mean_motion= SatList[i].meanmo #rev/day
    # M0_mean Anomally is SatList[i].meanan in degrees
    # n_dot_mean_motion is ndot/2 from the TLE data
    # n_2dots_mean_motion is n2dot/6 taken directly form TLE data
    
    
    
    #Assume Reference Epcoh is in TLE format
    refepoch=str(ts_sat_epoch)
    
    #Converts Reference Epoch Time to a Datetime Objects
    Epochdt=refepoch_to_dt(refepoch)
    
    
    Epochdt_list.append(Epochdt)
    
    
    t=(time-Epochdt).total_seconds()
    #time since Epoch in TLE
    
    time_since_epoch_sec.append(t)
    
    
    
    
    Mt_mean_anomaly=M0_mean_anomaly+ \
        n_mean_motion*360*((t/86400))+360*(n_dot_mean_motion)*((t/86400)**2)+ \
        360*(n_2dots_mean_motion)*((t/86400)**3)
        #outputs in deg
        
    Nt_mean_anomaly_motion=n_mean_motion* \
        (360/86400) + 2*360*(n_dot_mean_motion)*(t/(86400**2))+ \
        3*360*(n_2dots_mean_motion)*(t**2/(86400**3))
        # outputs in degs/sec
        
    Nt_mean_anomaly_motion_rev_day=Nt_mean_anomaly_motion*240    

    #Removing Mutlples
    Mt_mean_anomaly=Mt_mean_anomaly%360
    
    # Low Level Debug Helper for testing Mean Motion 
    global zTest_Mt_Mean_anomaly, zTest_Nt_mean_anomaly_motion_rev_day
    
    try:
        
        zTest_Mt_Mean_anomaly.append(Mt_mean_anomaly)
        zTest_Nt_mean_anomaly_motion_rev_day.append(Nt_mean_anomaly_motion_rev_day)
    except:
        zTest_Mt_Mean_anomaly=[]
        zTest_Mt_Mean_anomaly.append(Mt_mean_anomaly)
        zTest_Nt_mean_anomaly_motion_rev_day=[]
        zTest_Nt_mean_anomaly_motion_rev_day.append(Nt_mean_anomaly_motion_rev_day)
        
    
    
    return Mt_mean_anomaly,Nt_mean_anomaly_motion_rev_day

    # In[]
def KeplerEqn(Mt_mean_anomaly,eccentricity):
    
    #Examples Permitted Error
    permitted_error=0.05
    #Example Permitted Error
    
    
    Mt_mean_anomaly=np.deg2rad(Mt_mean_anomaly)
    #Converts to Radians
    
    e=float(eccentricity)
    #ensures that e is in the right format
    
    
    
    #Initialize lists
    E_=[]
    Del_M_=[]
    Del_E_=[]
    i=0
    
    #Calculates First Iteration
    E_.append(float(Mt_mean_anomaly))
    Del_M_.append(float(E_[i])-e*math.sin((E_[i]))-float(Mt_mean_anomaly))
    Del_E_.append(Del_M_[i]/(1-e*math.cos(E_[i])))
    Del_E_mag=abs(Del_E_[i])
    E_.append(E_[i]+Del_E_[i])
    
    #Calculates Further Iterations
    while Del_E_mag > permitted_error:
        
        Del_M_.append(float(E_[i])-e*math.sin((E_[i]))-float(Mt_mean_anomaly))
        Del_E_.append(Del_M_[i]/(1-e*math.cos(E_[i])))
        Del_E_mag=abs(Del_E_[i])
        (Del_E_mag)
        E_.append(E_[i]+Del_E_[i])
        i=i+1
        

   
    
        
    return np.degrees(E_[i+1]%(2*math.pi)) #reduces Eccentric Anom 
#! Returns in degrees
    # In[]
def perifocal(eccentricity,ecc_anomaly,a_semi_major_axis,omega_longitude_ascending_node, \
              omega_argument_periapsis,inclination,nt_mean_motion):
    
    #Assuming Earth
    mu=398600.4418 #km^3/s^2
    
    #Ensuring Orbital ELements are in right format
    eccentricity=float(eccentricity) 
    ecc_anomaly=float(ecc_anomaly)
    omega_longitude_ascending_node=float(omega_longitude_ascending_node)# deg
    omega_argument_periapsis=float(omega_argument_periapsis) #deg
    inclination=float(inclination) #deg
    nt_mean_motion=float(nt_mean_motion) #deg/sec
    
    ecc_anomaly=np.deg2rad(ecc_anomaly)
    
    
    #Calculating True Anomaly
    true_anom=2*(np.arctan(math.sqrt((1+eccentricity)/ \
                                     (1-eccentricity)))*np.tan(ecc_anomaly/2))
    true_anom=((np.arctan2(math.sin(true_anom),math.cos(true_anom)))+2*math.pi)%(2*math.pi)
    
        #radians
    # different approach
    v=math.acos((math.cos(ecc_anomaly)-eccentricity)/(1-eccentricity*math.cos(ecc_anomaly)))
        
      
        
        
    #Calculating R and its components
    r=a_semi_major_axis*(1-eccentricity**2)/(1+eccentricity*math.cos(true_anom))
    r_px=r*math.cos(true_anom)
    r_py=r*math.sin(true_anom)
    r_pz=0
    R_per=[r_px,r_py,r_pz]
    #Calculating Velocity Components
        #Note: We could not differentiate r in python so we used a work around method
        
    semi_lactus_rectum=a_semi_major_axis*(1-eccentricity**2)
    angular_mo=math.sqrt(mu*semi_lactus_rectum)
    v_px=-mu*math.sin(true_anom)/angular_mo
    v_py=mu*(eccentricity+math.cos(true_anom))/angular_mo
    v_pz=0
    v_per=[v_px,v_py,v_pz]
    #km/s
    
    # Low Level Debug Helper
    global zTest_true_anom,zTest_R_per, zTest_v_per
    try:
        zTest_R_per.append(R_per)
        zTest_v_per.append(v_per)
        zTest_true_anom.append(np.degrees(true_anom))
    except:
        zTest_R_per=[]
        zTest_v_per=[]
        zTest_true_anom=[]
        zTest_R_per.append(R_per)
        zTest_v_per.append(v_per) 
        zTest_true_anom.append(np.degrees(true_anom))
        
    
    return R_per,v_per
    # In[]
def sat_ECI(eccentricity,ecc_anomaly,a_semi_major_axis,omega_longitude_ascending_node,omega_argument_periapsis,inclination,nt_mean_motion):
    #Perifocal to ECI
    #Earth Centred Inertial Frame
    
    #Finds Perifocal Components
    r_per,v_per=perifocal(eccentricity,ecc_anomaly,a_semi_major_axis,omega_longitude_ascending_node,omega_argument_periapsis,inclination,nt_mean_motion)
    
    
    #Creates transformation
    Per_to_ECI=R.from_euler('ZXZ',[float(omega_longitude_ascending_node),float(inclination),float(omega_argument_periapsis)],degrees=True).as_matrix()
    #Note: RAAN,omega,inc in degrees
    pos_ECI = (np.matmul(Per_to_ECI,r_per)).tolist()
    
    #pos_ECI=(Per_to_ECI.apply(r_per)).tolist()
    vel_ECI=(np.matmul(Per_to_ECI,v_per)).tolist()
    #List variables are easier for me to deal with 
    
    
    
    
    
    return pos_ECI,vel_ECI
    # In[]
def sat_ECF(theta_t,eci_position,eci_velocity):
    
    
    #Creates rotational transformation
    ECI_to_ECF=R.from_euler('Z',[-float(theta_t)], degrees=True).as_matrix()
    
    #Applies Rotational Transformation
    #pos_ECF=ECI_to_ECF.apply(eci_position).tolist()[0] # had to perform weird
    
    pos_ECF=(np.matmul(ECI_to_ECF,eci_position)).tolist()[0]
    #weird conversion to get back to list
    #vel_ECF=ECI_to_ECF.apply(eci_velocity).tolist()[0]
    vel_ECF=(np.matmul(ECI_to_ECF,eci_velocity)).tolist()[0]
    #km/s
    
    
    #Relative Velocity
    Siderial_rotation=np.array([[0,0,0],[0,0,0],[0,0,360/86164.091]]) #Degrees/s
    
    #vel_rel=ECI_to_ECF.apply(eci_velocity-np.matmul(Siderial_rotation,eci_position))
    vel_rel=np.matmul(ECI_to_ECF,eci_velocity-np.matmul(Siderial_rotation,eci_position)).tolist()
    #km/s
    
    
    
    
    return pos_ECF,vel_ECF,vel_rel
    # In[]
def station_ECF(station_longitude,station_latitude,station_elevation):
    #input as rads only
    
    f=1/298.25223563
    R_e=6378.137 #km
    #geodetic longitde must be in degs
    #geodetic latitude must be in degs
    #station elevation must be in km
    phi=(float(station_latitude)*math.pi/180)
    h=float(station_elevation)
    lambda_=float(station_longitude)*math.pi/180
    e=math.sqrt(2*f-f**2)
    n_phi=R_e/(math.sqrt(1-(e**2)*(math.sin(phi))**2))
    
    
    T_x=(n_phi+h)*math.cos(phi)*math.cos(lambda_)
    T_y=(n_phi+h)*math.cos(phi)*math.sin(lambda_)
    T_z=(((1-e**2)*n_phi+h)*math.sin(phi))
    
    #R_x=station_body_position[0]-T_x
    #R_y=station_body_position[1]-T_y
    #R_z=station_body_position[2]-T_z
    
    #where station_body_poisitonh is the sat_ECF coordinates
    
    
    
    
    return [T_x,T_y,T_z]
    
    
  #Note: Professor said that we did not have to do the second station_EFC fucntion
    # In[]


def range_ECF2topo(station_body_position, \
                   sat_ecf_position,sat_ecf_velocity,station_longitude, \
                   station_latitude):
    
    station_longitude=np.deg2rad(float(station_longitude))
    station_latitude=np.deg2rad(float(station_latitude))
    
    
    
    
    #assuming that station body positon is [Tx,Ty,Tz]
    
    R=[sat_ecf_position[0]-station_body_position[0],\
       sat_ecf_position[1]-station_body_position[1],\
        sat_ecf_position[2]-station_body_position[2]]
    
    
    #Intializes Transformation Matrix
    T_ECF_to_topo=[[-math.sin(station_longitude), \
                    math.cos(station_latitude),0], \
                   [math.cos(station_longitude)*math.sin(station_latitude), \
                    -math.sin(station_longitude)*math.sin(station_latitude), \
                    math.cos(station_latitude)],[math.cos(station_longitude) \
                                                 *math.cos(station_latitude), \
                                                 math.sin(station_longitude)* \
                                                 math.cos(station_latitude), \
                                                 math.sin(station_latitude)]]
    
    
    #Transform Range Vector
    R_ti=np.matmul(np.array(R),np.array(T_ECF_to_topo))
    
    # Assuming that sat_ecf Velocity is Relative
    vel_rel=sat_ecf_velocity
    
    #Transform Velocity Vector
    v_rel_ti=np.matmul(np.array(vel_rel),np.array(T_ECF_to_topo)).tolist()[0]
    
    return R_ti,v_rel_ti
    # In[]
def range_topo2look_angle(range_topo_position,range_topo_velocity):
    R=range_topo_position
    #Assuming range_topo_velocity is relative
    v_rel=range_topo_velocity
    
    #Calculates the AZ and EL
    AZ=(np.arctan2(R[0],R[1])+2*math.pi)%(2*math.pi)
    EL=math.atan(R[2]/(math.sqrt(R[0]**2+R[1]**2))) #Range is 0->90
    
    
    r=np.linalg.norm(R) #scalar of R
    
    R_xy=[R[0],R[1]]
    #In Software Specification Rxy is [tx ty]{Rtx;Rty} which would give a result of a singular value
    #Here we assume the Professor meant R_xy= the x and y components of R
    v_xy=[v_rel[0],v_rel[1]]
    
    #Calculates rates of AZ and EL
    rate_of_AZ=np.cross(v_xy,R_xy)/((np.linalg.norm(R_xy))**2)
    rate_of_EL=(r*v_rel[2]-R[2]*np.dot(R_xy,v_xy)/r)/(r**2)
    
    return AZ,EL,rate_of_AZ,rate_of_EL

    # In[]


 ##                         Initialization Codes
     # In[] 

class Station():
    def __init__(self,STNFIL):
        STNfile=open(STNFIL,'rt')
        #opens file using directory extension
        self.name=STNfile.readline()
        self.stnlat=STNfile.readline()
        self.stnlong=STNfile.readline()
        self.stnalt=STNfile.readline()
        self.utc_offset=STNfile.readline()
        self.az_el_nlim=STNfile.readline()
        self.az_el_lim=[]
        Iteration=int(self.az_el_nlim)
        i=0
        while i<Iteration:
            self.az_el_lim.append(STNfile.readline())
            i=i+1
            
        #self.az_el_lim=STNfile.readline()
        self.st_az_speed_max=STNfile.readline()
        self.st_el_speed_max=STNfile.readline()
        #reads the file line by line
        
        STNfile.close
        #closes files

    # In[]
class Satellite():
    def __init__(self,line0,line1,line2):
        self.name=line0
        self.refepoch=line1[18:32]
        self.incl=line2[8:16]
        self.raan=line2[17:25]
        self.eccn="."+line2[26:33]
        self.argper=line2[34:42] #degrees
        self.meanan=line2[43:51] #degrees
        self.meanmo=line2[52:63] #rev per day 
        self.ndot=line1[33:43] #is actually nodot/2
        self.n2dot=line1[44:50] #is actually n2dot/6
        self.bstar=line1[53:61]
        self.orbitnum=line2[63:68]
    # In[]
#Propagates the SatList with Satellites from TLE file
def SatListPropagate(SatFIL):
    Satfile=open(SatFIL,'rt')
    entries=(len(open(SatFIL).readlines()))/3
    i=0
    while i<entries:
        line0=Satfile.readline()
        line1=Satfile.readline()
        line2=Satfile.readline()
        try:
            SatList.append(Satellite(line0,line1,line2))
            print(SatList[i].name,'Satellite has been registered as Satellite:',len(SatList),"Array index: [",len(SatList)-1,"]")
        except:
            SatList=[]
            SatList.append(Satellite(line0,line1,line2))
            print('SatList has been created with',SatList[i].name,"Being the First Satellite")
        i=i+1
    return SatList
    # In[]
class tracking():
    def __init__(self, TrackingData):
        Trackingdatafile=open(TrackingData,'rt')
        self.starttime=(Trackingdatafile.readline())[0:19]
        self.endtime=(Trackingdatafile.readline())[0:19]
        self.timestep=Trackingdatafile.readline()
        #It is assumed that the step time is in seconds
        Trackingdatafile.close   
        # In[]
class linkinput():
    def __init__(self, LinkInputs):
        Linkinputfile=open(LinkInputs, 'rt')
        self.frequency=Linkinputfile.readline()
        self.Antennaeff=Linkinputfile.readline()
        self.AntennaDia=Linkinputfile.readline()
        self.Bandwidth=Linkinputfile.readline()
        self.RCVgain=Linkinputfile.readline()
        self.RCVnoise=Linkinputfile.readline()
        Linkinputfile.close
        #t=4*math.pi*(46)/(3.0e8*1.575e9)
        #signalloss=20*math.log(10,4*math.pi*(46)/(3.0e8*1.57542e9))
        #Will print out signal loss from the given values. 46 is antenna diameter
        #3.0e8 is speed of light and 1.57542e9 is frequency in Hz from MHz
        #print(signalloss)
        

    # In[]

##                  Main Calling Functions    
        
        
        
        
    # In[]

#This is a high-level function which gathers the string associated with the directories 
# and access lower level functions to create Instances of each 
def User_Input_parser_Call(StationLocationStr,TLEfilestr,TrackingSchedulestr,LinkInputsstr):
    StationInstance=Station(StationLocationStr)
    #Creates a Station Instance by using the Station file directory
    
    SatList=SatListPropagate(TLEfilestr)
    #Creates a List of Satellites 
    
    Tracking=tracking(TrackingSchedulestr)
    LinkData=linkinput(LinkInputsstr)
    return StationInstance,SatList,Tracking,LinkData
    # In[]
def Sat_pos_velCall(StationInstance,SatList,Tracking):
    #StationInstant, Times, and Links are class instances with their own attributes
    # Inputs: Station Instance Calculated, Satellite List, Time 
    
    #Low Level Debug Helper
    global zTest_Ecc_anom,zTest_ECI_R,zTest_ECI_v,zTest_ECF_R,zTest_ECF_vel
    global zTest_ECF_vel_rel,zTest_T_x,zTest_T_y,zTest_T_z
    
    zTest_Ecc_anom=[]
    zTest_ECI_R=[]
    zTest_ECI_v=[]
    zTest_ECF_R=[]
    zTest_ECF_vel=[]
    zTest_ECF_vel_rel=[]
    zTest_T_x=[]
    zTest_T_y=[]
    zTest_T_z=[]
    
    #We are assuming that in the next version of the code we will be iterating through Start and End times 
    #With a time step. This will be our Time value
    Time_start_dt=dt.datetime.strptime(Tracking.starttime,'%Y-%m-%d-%H:%M:%S') #This is to be used as a place holder until we iterate through the times 
    Time_end_dt=dt.datetime.strptime(Tracking.endtime,'%Y-%m-%d-%H:%M:%S')
    
    #Assuming that timesteps is in seconds 
    Time_iterations=(Time_end_dt-Time_start_dt).total_seconds()/float(Tracking.timestep)
    
    
    Time_dt=Time_start_dt
    #initializing empty arrays
    time=[]
    Satname=[]
    AZ_list=[]
    EL_list=[]
    Rate_of_AZ_list=[]
    Rate_of_EL_list=[]
    R_ti_list=[]
    v_rel_ti_list=[]
    Satnum_list=[]
    #Satnum List helps with identifying the Satellite
    
    global Epochdt_list
    Epochdt_list=[]
    global time_since_epoch_sec
    time_since_epoch_sec=[]
    global a_list
    a_list=[]
    
    #Propagtes data for THETAN
    Refepoch=referenceepoch_propagate(Tracking)
        #THETAN has been edited to input time_start_dt and Time_dt
    GMST=THETAN(Refepoch)
    
    global zTest_GMST_List
    zTest_GMST_List=GMST

#Iterates through satellite list first then for time 
#creates list of 
    for i in range(0,int(Time_iterations)):
      
      for p in range(0,int(len(SatList))):
        
        
        [Mt_Mean_anomaly,Nt_anomaly_motion]=mean_anomaly_motion(Time_dt,SatList[p].refepoch,float(SatList[p].meanan),float(SatList[p].meanmo),float(SatList[p].ndot),float(SatList[p].n2dot))
        #degrees,rev/day
        ecc_anomaly=KeplerEqn(Mt_Mean_anomaly,SatList[p].eccn)
        
        
        
        mu=398600.4418 #km^3/s
        
        a=(mu/(2*np.pi*float(Nt_anomaly_motion)/86400)**2)**(1/3)
        a_list.append(a)
        
        [pos_ECI,vel_ECI]=sat_ECI(SatList[p].eccn,ecc_anomaly, \
        a,SatList[p].raan,SatList[p].argper,SatList[p].incl,Nt_anomaly_motion)
        
        
        GMST_1=GMST[i]    
        [pos_ECF,vel_ECF,vel_rel_ECF]=sat_ECF(GMST_1,pos_ECI,vel_ECI)
        #Note: We assume that station_body_position is Tx,Ty,Tz
        [Tx,Ty,Tz]=station_ECF(StationInstance.stnlong,StationInstance.stnlat,StationInstance.stnalt)
    
    
         #Note: Station Long and Latitude must be in Radians
        [R_ti,v_rel_ti]=range_ECF2topo([Tx,Ty,Tz],pos_ECF,vel_rel_ECF,StationInstance.stnlong,StationInstance.stnlat)
        
        
        [AZ,EL,Rate_of_AZ,Rate_of_EL]=range_topo2look_angle(R_ti,v_rel_ti)



        Satnum_list.append(p)
        AZ_list.append(AZ)
        EL_list.append(EL)
        Rate_of_AZ_list.append(Rate_of_AZ)
        Rate_of_EL_list.append(Rate_of_EL)
        R_ti_list.append(R_ti)
        v_rel_ti_list.append(v_rel_ti)

        time.append(Time_dt)
        
        
        # Low Level Debugging Helper
        zTest_Ecc_anom.append(ecc_anomaly)
        
        zTest_ECI_R.append(pos_ECI),zTest_ECI_v.append(vel_ECI)
        zTest_ECF_R.append(pos_ECF),zTest_ECF_vel.append(vel_ECF),zTest_ECF_vel_rel.append(vel_rel_ECF)
        zTest_T_x.append(Tx),zTest_T_y.append(Ty),zTest_T_z.append(Tz)
        
        
        
        
        
      Time_dt=Time_dt+dt.timedelta(seconds=float(Tracking.timestep))
      #At the End of the loop add time 
      
    
    
    #Start Time is in EST and is converting inside function
    #Note: we have made changes to the THETAN code to also input Tracking.starttime. As of this version, 
    #Time now is used for the t variable. This is to be changed in later versions when we iterate through time 
    
    
    
    
    
    
    
    
    return AZ_list,EL_list,Rate_of_AZ_list,Rate_of_EL_list,R_ti_list,v_rel_ti_list,time,Satnum_list

    # In[]
def Pointing(StnInstance,AZ_list,EL_list,time,Satnum_list,Signal_lost):
  
  
  #avail=available for viewing
  AZ_avail=[]
  EL_avail=[]
  Times_avail=[]
  Satnum_avail=[]
  #creates blanks arrays to be filled
  AZ_AOS=[]
  EL_AOS=[]
  Times_AOS=[]
  SatNum_AOS=[]
  AZ_LOS=[]
  EL_LOS=[]
  Times_LOS=[]
  SatNum_LOS=[]
  Satnum_iteration=max(Satnum_list)+1
  time_delta=time[1]-time[0]
  Signal_lost_AOS=[]
  Signal_lost_LOS=[]
  Signal_lost_avail=[]
  
  #this creates a variable that changes when Available, AOS, or LOS
  #it's the same length as AZ, EL, times to benefit CSV file
  global Avail_list, AOS_List_boolean, LOS_List_boolean
  Avail_list=[]
  AOS_List_boolean=[]
  LOS_List_boolean=[]
  for q in range(0,len(AZ_list)):
      Avail_list.append(0)
      AOS_List_boolean.append(0)
      LOS_List_boolean.append(0)
      
  #
  
  

  
  for i in range(0,int(StnInstance.az_el_nlim)):
      
      # For Each iteration of the Station Instation Limits
    ThisStationLimit=StnInstance.az_el_lim[i].split(",")
    
    # Covnerts Limits from Degrees to Rads
    AZ_muth_max=np.deg2rad(float(ThisStationLimit[0]))
    if i>0:
        PreviousStationLimit=StnInstance.az_el_lim[i-1].split(",")
        AZ_muth_min=np.deg2rad(float(PreviousStationLimit[0]))
        #min Value for this Azimuth Step
    else:
        AZ_muth_min=0
        if(int(StnInstance.az_el_nlim)==1):
            AZ_muth_max=math.pi*2
        
    
    EL_lim_max=np.deg2rad(float(ThisStationLimit[2]))
    EL_lim_min=np.deg2rad(float(ThisStationLimit[1]))
    
    
    for j in range(0,len(AZ_list)):
      
      if AZ_list[j] > (AZ_muth_min) and AZ_muth_max > AZ_list[j] and (EL_lim_max) > EL_list[j] and EL_list[j] > (EL_lim_min):
        #compares Azimuth and Elevations to Limits
          
          #Satellite Available
          
        
          #Compares to Previous Value
          # If previous value is not in limits but is now either AOS or from another limit iteration** 
        if j >= (Satnum_iteration) and  Satnum_list[j] == Satnum_list[j-Satnum_iteration] and (AZ_list[j-Satnum_iteration] < float(AZ_muth_min) or AZ_list[j-Satnum_iteration] > float(AZ_muth_max) or float(EL_lim_max) < EL_list[j-Satnum_iteration] or EL_list[j-Satnum_iteration] < float(EL_lim_min)):
            # A Signal has been acquired
            
            #** Here we test to see if the AOS was caused by limit bound issues or from actual AOS
            
            if i > 0:
                
                    
                
                
                
                append=0  
                for k in range(0,len(Times_LOS)):
                  
                    
                    if time[j] == Times_LOS[k] and Satnum[j]== SatNum_LOS[k]:
                        append=0
                        #False Signal Acquire- Delete previous LOS
                        del Times_LOS[k]
                        del SatNum_LOS[k]
                        del AZ_LOS[k]
                        del EL_LOS[k]
                        LOS_List_boolean[k]=0
                    else: 
                        append=1
                        #True Signal Acquired
                        
                        break
                if append==1:
                    AZ_AOS.append(AZ_list[j])
                    EL_AOS.append(EL_list[j])
                    Times_AOS.append(time[j])
                    SatNum_AOS.append(Satnum_list[j])
                    Signal_lost_AOS.append(Signal_lost[j])
                    AOS_List_boolean[j]=1
                
            else:
                # First Signal Acquisition of Signal    
                AZ_AOS.append(AZ_list[j])
                EL_AOS.append(EL_list[j])
                Times_AOS.append(time[j])
                SatNum_AOS.append(Satnum_list[j])
                Signal_lost_AOS.append(Signal_lost[j])
                
                
                
                 
                        
                        
                
                
             
                
                
                
                
        AZ_avail.append(AZ_list[j])
        EL_avail.append(EL_list[j])
        Times_avail.append(time[j])
        Satnum_avail.append(Satnum_list[j])
        Signal_lost_avail.append(Signal_lost[j])
        Avail_list[j]=1
      else:
          # Satellite Unavailable
          
          #Compares previous iteration of Satellite to current Satellite Available
          if j >= (Satnum_iteration)  and Satnum_list[j]==Satnum_list[j-(Satnum_iteration)]  and AZ_list[j-(Satnum_iteration)] > (AZ_muth_min) and AZ_muth_max > AZ_list[j-(Satnum_iteration)]  and float(EL_lim_max) > EL_list[j-(Satnum_iteration)] and EL_list[j-(Satnum_iteration)] > float(EL_lim_min):
              #If Satellite was available but is no longer then it is either LOS or out of our Limit Range
              #
              
              AZ_LOS.append(AZ_list[j])
              EL_LOS.append(EL_list[j])
              Times_LOS.append(time[j])
              SatNum_LOS.append(Satnum_list[j])
              Signal_lost_LOS.append(Signal_lost[j])
              LOS_List_boolean[j]=1
              
              
              

    
  AOS_List=[AZ_AOS,EL_AOS,Times_AOS,SatNum_AOS,Signal_lost_AOS]
  LOS_List=[AZ_LOS,EL_LOS,Times_LOS,SatNum_LOS,Signal_lost_LOS]
#Creates a list of Available Azimuth, Elevation, Times and Satnum available 
  
  
  return AZ_avail,EL_avail,Times_avail,Satnum_avail,AOS_List,LOS_List


    # In[]
#def Link_Calculations(LinkData):


def linkInputscall(linkdat):
  Linkcallfile=open(linkdat, 'rt')
  frequency=Linkcallfile.readline()
  Antennaeff=Linkcallfile.readline()
  AntennaDia=Linkcallfile.readline()
  Linkcallfile.close
  #signalloss=20*math.log(10,4*math.pi*(float(AntennaDia))/(3.0e8*float(frequency)*1e6))
  return frequency,Antennaeff,AntennaDia

    # In[]
def Visibility(StationInstance,AZ,EL,times,Satnum,Signal_lost):
    #[AZ_avail,EL_avail,Times_avail,Satnum_avail]=Pointing(StationInstance,AZ,EL,Satnum)
    AOS=Pointing(StationInstance,AZ,EL,times,Satnum,Signal_lost)[4]
    LOS=Pointing(StationInstance,AZ,EL,times,Satnum,Signal_lost)[5]
    Satnum_AOS=AOS[3]
    Satnum_LOS=LOS[3]
    Sat_Signal_Lost=AOS[4]
    
    Sat_AOS_Time=AOS[2]
    Sat_LOS_Time=LOS[2]
    
    AOS_LOS_list=[]
    
    
            
    
    
    
    for i in range(0,len(Satnum_AOS)):
        for j in range(0,len(Satnum_LOS)):
            
            if Satnum_AOS[i] == Satnum_LOS[j] and Sat_LOS_Time[j] >= Sat_AOS_Time[i] :
                
                Templist=[Satnum_AOS[i],SatList[Satnum_AOS[i]].name,Sat_AOS_Time[i],Sat_LOS_Time[j],Sat_Signal_Lost[i]]
                #creates a temporary list to store Sat number,Sat list name ,AOS time,LOS time and 
                AOS_LOS_list.append(Templist)
    

        
    global Unique_AOS_LOS
    Unique_AOS_LOS=[[],[],[],[],[]]
    
   
    for i in range(0,len(AOS_LOS_list)):
        AOS_time=AOS_LOS_list[i][2]
        Satname=AOS_LOS_list[i][1]
        Satnum=AOS_LOS_list[i][0]
        LOS_time=AOS_LOS_list[i][3]
        Sat_Signal_Lost=AOS_LOS_list[i][4]
        if AOS_time not in Unique_AOS_LOS[2]:
            
            #creates a Unique List of AOS
            Unique_AOS_LOS[0].append(Satnum)
            Unique_AOS_LOS[1].append(Satname)
            Unique_AOS_LOS[2].append(AOS_time)
            Unique_AOS_LOS[3].append(LOS_time)
            Unique_AOS_LOS[4].append(Sat_Signal_Lost)
            
            
    
    
    return Unique_AOS_LOS
    # In[]
def SingalCalc(freq,Antennaeff,AntennaDia,R_ti):
    Signal_loss=[]
    DopplerShift_list=[]
    for i in range(0,len(R_ti)):
        R=np.linalg.norm(R_ti[i])
        Signal_loss.append(20*math.log(10,4*math.pi*(float(R))/(3.0e8*float(freq)*1e6)))
        v=np.dot(v_rel_ti[i],R_ti[i])/R
        #output in km
        c=(2.998e+8)/1000 #m/s->km/s
        DopplerShift_list.append(-v*float(freq)*1e6/c) #converts frequency to Hz
    return Signal_loss,DopplerShift_list

 
   

    # In[]
#       Debugging 
    
    
    # In[]
def STKout(EphemFile,Epochtime,time,Coord,position,velocity):
    #Writes Ephermis File
    
    
    #Asssuming StartString is the Start Time String
    
    #Here I interpreted the positon and velocity arrays as
    # Postion Array:
    #  X Y Z
    
    #opens or creates a file with the name EphemFile
    file_1=open(EphemFile,"w+")
    
    #Creates a set of Strings that can be easily changed. 
    
    ScenarioEpochString="ScenarioEpoch \t"+Epochtime
    #Scenario Epch is not requried 
    CentralBody="Earth"
    CentralBodyString="CentralBody "+ CentralBody
    CoordinateSystemString="CoordinateSystem "+Coord
    NumofEphermis=len(time) #I don't know how to attain this number
    
    #Writing Header
    file_1.write("stk.v.5.0 \nBEGIN Ephemeris \n")
    
    file_1.write("NumberOfEphemerisPoints "+str(NumofEphermis)+"\n")
    file_1.write(ScenarioEpochString+"\n")
    file_1.write("DistanceUnit \t Kilometers \n")
    file_1.write(CentralBodyString+"\n")
    file_1.write(CoordinateSystemString+"\n")
    file_1.write("EphemerisTimePosVel\n")
    
    #extracts values from position and velocity lists
    X=[row[0] for row in position]
    Y=[row[1] for row in position]
    Z=[row[2] for row in position]
    X_dot=[row[0] for row in velocity]
    Y_dot=[row[1] for row in velocity]
    Z_dot=[row[2] for row in velocity]
    
    i=0
    for i in range(0,len(time)):
        
        file_1.write('{0:15.14e} {1:15.14e} {2:15.14e} {3:15.14e} {4:15.14e} {5:15.14e} {6:15.14e} \n'.format(float(time[i]),float(X[i]),float(Y[i]),float(Z[i]),float(X_dot[i]),float(Y_dot[i]),float(Z_dot[i])))
        #writing to file in formatted way
        
        
        
        # for i is in time, we iterate through the list to write in the values to the value
        i=i+1
        
    file_1.write("\n\nEND Ephemeris")
    file_1.close()

# In[]
def STKsp(AZ,EL,time,SpFile):
   #This function outputs the Pointing File used with STK    
    
    #Takes in Aximuth and Elevatyion in Radians
        #We should have an input file name. I've provided one here to help
        
        
        #opens file
        file_2=open(SpFile,"w+")
        
        NumofAttitudes=len(time) 
        
        #Writing Header
        file_2.write("stk.v.4.3 \nBEGIN\tAttitude\n")
        file_2.write("NumberOfAttitudePoints "+str(NumofAttitudes)+"\n")
        file_2.write("AttitudeTimeAzElAngles\n")
        
        
    
        for i in range(0,len(time)):
            file_2.write('{0:7.2f} {1:7.2f} {2:7.2f}\n'.format(float(time[i]),float(np.rad2deg(AZ[i])),float(np.rad2deg(EL[i]))))
        # for i is in time, we iterate through the list to write in the values to the value
            
            
        #writing closer 
        file_2.write("END Attitude")
        file_2.close    
        
        # In[]
#This function write the Mater.csv file which is used extensively for Debugging
def Master_csvwriter(filename,AZ,EL,Rate_of_AZ,Rate_of_EL,Mean_anomaly,Mean_anomaly_motion,R_ti,v_rel_ti,time,Satnum,Avail_list,AOS_List_boolean,LOS_List_boolean):
    with open(filename,mode='w+',newline='') as csv_file:
        
        csv_writer=csv.writer(csv_file,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(["Perifocal Range","Perifocal Velocity","ECI Position","ECI_Velocity","ECF Position","ECF Velocity" \
                             ,"Azimuth","Elevation","Rate of AZ","Rate of EL","Mean Anom @time","Mean Anom Motion @time","Topocentric Range",\
                                 "Topocentric Reletive Velocity","Time","Time Since Epoch","Satellite Name","Avilable for View TRUE=1", \
                                     "Signal Acquired","Signal Lost"])
        for s in range(0,len(AZ)):
           csv_writer.writerow([zTest_R_per[s],zTest_v_per[s],zTest_ECI_R[s],zTest_ECI_v[s],zTest_ECF_R[s],zTest_ECF_vel_rel[s],AZ[s],EL[s],Rate_of_AZ[s],Rate_of_EL[s],Mean_anomaly[s],Mean_anomaly_motion[s],R_ti[s],v_rel_ti[s],time[s],time_since_epoch_sec[s],SatList[Satnum[s]].name,Avail_list[s],AOS_List_boolean[s],LOS_List_boolean[s]])
       
        csv_file.close    
        return 
# This functions converts all the lists into a csv file which can be read into excel and understood easily 


 # In[]
# Note there are some issues involving the way AOS_LOS_list is orientated. 
def AOS_csvwriter(filename,AOS_LOS_List):
    with open(filename,mode='w+',newline='') as csv_file:
        
        csv_writer=csv.writer(csv_file,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(["Sat No.","Sat Name","AOS Time","LOS Time","Min Expected Level",])
        for s in range(0,len(AOS_LOS_List[0])):
            #AOS_instance=AOS_LOS_List[s]
            csv_writer.writerow([AOS_LOS_List[0][s],AOS_LOS_List[1][s],AOS_LOS_List[2][s],AOS_LOS_List[3][s],AOS_LOS_List[4][s]])
       
        csv_file.close    
        return 
    
# In[]
def AZ_EL_csvwriter(filename,Satnum_avail,AZ_avail,EL_avail,Time_Avail):
    with open(filename,mode='w+',newline='') as csv_file:
        
        csv_writer=csv.writer(csv_file,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(["Sat Name","Azimuth(rads)","Elevation(rads)","Time"])
        for s in range(0,len(Satnum_avail)):
            
            csv_writer.writerow([SatList[Satnum_avail[s]].name,AZ_avail[s],EL_avail[s],Time_Avail[s]])
       
        csv_file.close    
        return
# In[]
   
    
    # In[]

def ChosenSat_Output_function(position,velocity,time,t_list,Epochdt_list,AZ,EL,AZ_rate,EL_Rate,R_ti,v_rel_ti,Sat_dB,DopplerShift,InertialEpemFileName,FixedEphemFileName,spFileName,ControlFileName):
    

    val=input("Enter the Satellite Index to Output in .sp and .e file:")
    print(SatList[int(val)].name, "Has been chosen")
    Satnum_iteration=len(SatList)
    
    sat_time=[]
    time_since_epoch_sec_sat=[]
    sat_position=[]
    sat_velocity=[]
    sat_position_ECF=[]
    sat_velocity_ECF=[]
    sat_AZ=[]
    sat_EL=[]
    sat_AZ_rate=[]
    sat_EL_rate=[]
    sat_Range=[]
    sat_v_rel_ti=[]
    sat_dB=[]
    sat_DS=[]
    
    
    for i in range(0,len(zTest_ECF_vel_rel)):
        zTest_ECF_vel_rel[i]=zTest_ECF_vel_rel[i][:][0]
    Num_of_iterations=len(time)/Satnum_iteration
    for i in range(0,int(Num_of_iterations)):
        #Iterates through time and timesince epoch for a specific Satellite
        
        #This is done because the time and T_list for a specific sat is seperated by the lenght of a Satlist
        
        
        #Here we gather oall the data corresponding to the selected Satellite
        sat_time.append(time[int(val)+i*Satnum_iteration])
        time_since_epoch_sec_sat.append(t_list[int(val)+i*Satnum_iteration])
        sat_position.append(position[int(val)+i*Satnum_iteration])
        sat_velocity.append(velocity[int(val)+i*Satnum_iteration])
        sat_position_ECF.append(zTest_ECF_R[int(val)+i*Satnum_iteration])
        sat_velocity_ECF.append(zTest_ECF_vel_rel[int(val)+i*Satnum_iteration])
        sat_AZ.append((np.rad2deg(AZ[int(val)+i*Satnum_iteration])))
        sat_EL.append((np.rad2deg(EL[int(val)+i*Satnum_iteration])))
        sat_AZ_rate.append((np.rad2deg(AZ_rate[int(val)+i*Satnum_iteration])))
        sat_EL_rate.append((np.rad2deg(EL_Rate[int(val)+i*Satnum_iteration])))
        sat_Range.append((np.linalg.norm(R_ti[int(val)+i*Satnum_iteration])))
        sat_v_rel_ti.append((np.linalg.norm(v_rel_ti[int(val)+i*Satnum_iteration])))
        sat_dB.append((Sat_dB[int(val)+i*Satnum_iteration]))
                                                
        sat_DS.append((DopplerShift[int(val)+i*Satnum_iteration]/1000))



# Outputting STK files
        
    #Change to right format    
    EpochTimeString=dt.datetime.strftime(Epochdt_list[int(val)],"%d %b %Y %H:%M:%S")  
    STKout(InertialEpemFileName,EpochTimeString,time_since_epoch_sec_sat,"Inertial",sat_position,sat_velocity)
    
    
         
    STKout(FixedEphemFileName,EpochTimeString,time_since_epoch_sec_sat,"Fixed",(sat_position_ECF),(sat_velocity_ECF))

    STKsp(sat_AZ,sat_EL,time_since_epoch_sec_sat,spFileName)
    
    
    print ('\n UTC\t AZ Deg EL Deg AZ-Vel (deg/sec) El-Vel (deg/sec Range (km)/sec Doppler KHz Level dBm\n')
#Outputting Tracking Data   
    for j in range(0,len(sat_time)):
        DOY_=doy(sat_time[j].year, sat_time[j].month, sat_time[j].day)
        print('{0:4.0f}{1:1}{2:3}{3:1}{4:8} {5:3.0f} {6:2.0f} {7:3.1f} {8:3.1f}  {9:2.0f} {10:2.0f} {11:3.1f} {12:3.1f}\n'\
                       .format(sat_time[i].year,"-",DOY_,"-",sat_time[i].ctime()[11:19],int(sat_AZ[i]),(sat_AZ[i]*60)%60,\
                           (sat_AZ[i]*3600)%60,sat_AZ_rate[i],int(sat_EL[i]),(sat_EL[i]*60)%60,(sat_EL[i]*3600)%60,sat_EL_rate[i]))
    
   # ,sat_AZ,sat_EL,sat_AZ_rate,sat_EL_rate,sat_Range,sat_v_rel_ti,sat_dB
    TrackingDataYN=input("Is this data acceptable?(Y/N)")
    if TrackingDataYN=='Y':
        print('...Printing Control Data...')
        
        
        file=open(ControlFileName,"w+")
        
        
        
        #Writing Header
        file.write("#ARO CONTROL FILE\n")
        file.write("#Station: ARO\n")
        file.write("#UTC DATE/TIME\t Azimuth and AZ_Velocity     Elevation and EL_Velocity\n")
        
        
    
        for i in range(0,len(sat_time)):
            DOY_=doy(sat_time[i].year, sat_time[i].month, sat_time[i].day)
            #Claculates Day of Year for time instant
            
            file.write('{0:4.0f}{1:1}{2:3}{3:1}{4:8} {5:3.0f} {6:2.0f} {7:3.1f} {8:3.1f}     {9:2.0f} {10:2.0f} {11:3.1f} {12:3.1f}\n'\
                       .format(sat_time[i].year,".",DOY_,".",sat_time[i].ctime()[11:19],int(sat_AZ[i]),(sat_AZ[i]*60)%60,\
                           (sat_AZ[i]*3600)%60,sat_AZ_rate[i],int(sat_EL[i]),(sat_EL[i]*60)%60,(sat_EL[i]*3600)%60,sat_EL_rate[i]))
        # for i is in time, we iterate through the list to write in the values to the value
            
            
        #writing closer 
        print('Done!')
        file.close 
    else:
        print("Okay Goodbye")
    return 

 # In[]

##                  Calling Functions 
# These functions largly control the inputs and outputs of the program. 
    #Edit Directories for control over the proper inputs and outputs of the program

#Assuming That The Use has already initialized all necessary functions and classes
# The Main Program can be deduced to this
[StationInstance,SatList,Tracking,LinkData]=User_Input_parser_Call(r'D:\School\5th Year Fall Semester\ESSE 4350\Tracking\P6\Station.txt',r'D:\School\5th Year Fall Semester\ESSE 4350\Tracking\P6\gps-ops.txt',r'D:\School\5th Year Fall Semester\ESSE 4350\Tracking\P6\TrackingData.txt',r'D:\School\5th Year Fall Semester\ESSE 4350\Tracking\P6\LinkInputs.txt')
# This call function creates instances of each file that can be easily manipulated in Python.

#This function calls upon the Satellite position velocity calculator functions
[AZ,EL,Rate_of_AZ,Rate_of_EL,R_ti,v_rel_ti,time,Satnum]=Sat_pos_velCall(StationInstance,SatList,Tracking)


#Gathers Link Inputs and assigns them as variables
[freq,Antennaeff,AntennaDia]=linkInputscall(r'D:\School\5th Year Fall Semester\ESSE 4350\Tracking\P6\LinkInputs.txt')
#MHZ


#Calculates Signal loss
[Signal_loss,DopplerShift]=SingalCalc(freq,Antennaeff,AntennaDia,R_ti)
#





[AZ_avail,EL_avail,Times_avail,Satnum_avail,AOS_List,LOS_List]=Pointing(StationInstance,AZ,EL,time,Satnum,Signal_loss)

#Visibility creates a formatted list 
AOS_LOS_list=Visibility(StationInstance,AZ,EL,time,Satnum,Signal_loss)


#Outputs CSV Files
os.makedirs(os.path.dirname("D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/Master.csv"),exist_ok=True)

#writing to a csv file
Master_csvwriter("D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/Master.csv",AZ,EL,Rate_of_AZ,Rate_of_EL,zTest_Mt_Mean_anomaly,zTest_Nt_mean_anomaly_motion_rev_day,R_ti,v_rel_ti,time,Satnum,Avail_list,AOS_List_boolean,LOS_List_boolean)
AOS_csvwriter("D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/AOS_LOS.csv",AOS_LOS_list)
AZ_EL_csvwriter("D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/AZ_EL.csv",Satnum_avail,AZ_avail,EL_avail,Times_avail)
#Outputs AZ in Rads


# Outputs STK files
ChosenSat_Output_function(zTest_ECI_R,zTest_ECI_v,time,time_since_epoch_sec,Epochdt_list,AZ,EL,Rate_of_AZ,Rate_of_EL,R_ti,v_rel_ti,Signal_loss,DopplerShift,\
                    'D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/EphemFileExampleInertial.e',\
                        'D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/EphemFileExampleFixed.e',\
                            "D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/STKSP.sp",\
                                "D:/School/5th Year Fall Semester/ESSE 4350/Tracking/P6+/OutputFiles/ControlFile.txt")
    


    



