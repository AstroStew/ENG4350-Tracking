# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 18:54:39 2020

@author: mstew
"""
##                  Draft Version 1


 ##                                 Import Sections
import datetime as dt
import numpy as np
import math
from scipy.spatial.transform import Rotation as R



 ##                   Satellite Position,Velocity Calculator Functions



    # In[]

#Converts time in TLE format to Datetime Format
def refepoch_to_dt(refepoch):
    Epochyrday = dt.datetime.strptime((refepoch[0:5]),'%y%j')
    dfrac = np.modf(np.float(refepoch))[0]
    dfracdt = dt.timedelta(microseconds=np.int(dfrac*24*3600*10**6))
    Epochdt = Epochyrday + dfracdt
    return Epochdt

    # In[]
# Finds the Day of Year from the Year Month and Day parameters
def doy(YR,MO,D):
    if len(str(MO))==1:
        MO='0'+ str(MO)
    if len(str(D))==1:
        D='0'+ str(D)
    String=str(YR)+str(MO)+str(D)
    Time = dt.datetime.strptime(String,'%Y%m%d')
    #Converts First to dt
    DOY=dt.datetime.strftime(Time,'%j')
    #Converts from dt to day of the year
    return DOY
    
    # In[]
# This creates a list of 
def referenceepoch_propagate(TrackingData):
    starttime=dt.datetime.strptime(TrackingData.starttime,'%Y-%m-%d-%H:%M:%S')
    endtime=dt.datetime.strptime(TrackingData.endtime,'%Y-%m-%d-%H:%M:%S')
    Iterations=(endtime-starttime).total_seconds()/float(TrackingData.timestep)
    refepoch_list=[]
    time=starttime
    
    for i in range(0,int(Iterations)):
        yearstring=str(time.year)[2:4]
        dayfrac=round(((time.microsecond/8.6e+10)+(time.second/86400)+ \
                          (time.minute/1440)+(time.hour/24)),4)
        doy_nofrac=str(doy(time.year,time.month,time.day))
        doystring=str(int(doy(time.year,time.month,time.day))+dayfrac)[0:13]
        time=time+dt.timedelta(seconds=float(TrackingData.timestep))
        
        
        if len(doy_nofrac)==3:
            resultstring=yearstring+doystring
        elif len(doy_nofrac)==2:
            resultstring=yearstring+"0"+doystring
        elif len(doy_nofrac)==1:
            resultstring=yearstring+"00"+doystring
            
        
            
        
        
        
        refepoch_list.append(resultstring)
        
        
    return refepoch_list
    
    # In[]
def THETAN(refepoch):
    #Input is a refepoch array this dsoesn't make sense as Tracking Data contains an easily parsible 
    
    
    GMST_list=[]
    
    
    J2000=dt.datetime.strptime('2000-01-01 12:00:00','%Y-%m-%d %H:%M:%S')
    global Starttime
    Starttime=dt.datetime.strptime(Tracking.starttime,'%Y-%m-%d-%H:%M:%S')
    Midtime=Starttime.replace(hour=0,minute=0,second=0)
    D_u=(Midtime-J2000).days+(Midtime-J2000).seconds/86400
    T_u=D_u/36525
    GMST_00=(99.9677947+36000.7700631*T_u+0.00038793*T_u**2-2.6e-8*T_u**3)%360

    
    

    for i in range(0,len(refepoch)):
        times=(refepoch_to_dt(refepoch[i]))
        
        #Creates T mid for Observation Day
        #Notice how we replace hour,min and sec to 0. This makes the time midnight!
        del_sec=(times-Midtime).total_seconds()
        

        D_u_2=(times-J2000).days+(times-J2000).seconds/86400

        T_u_2=D_u_2/36525

        r=1.002737909350795+5.9006e-11*(T_u_2-T_u)-5.9e-15*(T_u_2-T_u)**2
        GMST_t=(GMST_00+360*r/86400*(del_sec))%360
        
        GMST_list.append(GMST_t)
    
    #Low Level Debug Helper
    
        
    return GMST_list    

    # In[]
def mean_anomaly_motion(time,ts_sat_epoch,M0_mean_anomaly,n_mean_motion, \
                        n_dot_mean_motion,n_2dots_mean_motion):
    
    #Assume Reference Epcoh is in TLE format
    refepoch=ts_sat_epoch
    Epochyrday = dt.datetime.strptime(refepoch[:4],'%y%j')
    dfrac = np.modf(np.float(refepoch))[0]
    dfracdt = dt.timedelta(microseconds=np.int(dfrac*24*3600*10**6))
    Epochdt = Epochyrday + dfracdt
    
    #assume Time is datetime object
    
    t=(time-Epochdt).total_seconds()
    
    
    
    Mt_mean_anomaly=M0_mean_anomaly+ \
        n_mean_motion*(360*t/86400)+360*(n_dot_mean_motion)*(t/86400)**2+ \
        360*(n_2dots_mean_motion)*(t/86400)**3
    Nt_mean_anomaly_motion=n_mean_motion* \
        (360*t/86400) + 2*360*(n_dot_mean_motion)*(t/86400**2)+ \
        3*360*(n_2dots_mean_motion)*(t**2/86400**3)

    #Removing Mutlples
    Mt_mean_anomaly=Mt_mean_anomaly%360
    
    # Low Level Debug Helper
    global zTest_Mt_Mean_anomaly
    
    try:
        
        zTest_Mt_Mean_anomaly.append(Mt_mean_anomaly)
    except:
        zTest_Mt_Mean_anomaly=[]
        zTest_Mt_Mean_anomaly.append(Mt_mean_anomaly)
        
        
    
    
    return Mt_mean_anomaly,Nt_mean_anomaly_motion

    # In[]
def KeplerEqn(Mt_mean_anomaly,eccentricity):
    
    #Examples Permitted Error
    permitted_error=0.32
    #Example Permitted Error
    
    
    Mt_mean_anomaly=float(Mt_mean_anomaly)*math.pi/180
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
        
     
   
    
        
    return E_[i+1]%(2*math.pi) #reduces Eccentric Anom 
#! Returns in Radians
    # In[]
def perifocal(eccentricity,ecc_anomaly,a_semi_major_axis,omega_longitude_ascending_node, \
              omega_argument_periapsis,inclination,nt_mean_motion):
    
    #Assuming Earth
    mu=398600.4418 #km^3/s^2
    
    #Ensuring Orbital ELements are in right format
    eccentricity=float(eccentricity) 
    ecc_anomaly=float(ecc_anomaly)
    omega_longitude_ascending_node=float(omega_longitude_ascending_node)
    omega_argument_periapsis=float(omega_argument_periapsis)
    inclination=float(inclination)
    nt_mean_motion=float(nt_mean_motion)
    
    #Calculating True Anomaly
    true_anom=2*(math.atan(math.sqrt((1+eccentricity)/ \
                                     (1-eccentricity))*math.tan(ecc_anomaly)))
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
    global zTest_R_per, zTest_v_per
    try:
        zTest_R_per.append(R_per)
        zTest_v_per.append(v_per)
    except:
        zTest_R_per=[]
        zTest_v_per=[]
        zTest_R_per.append(R_per)
        zTest_v_per.append(v_per) 
    
    return R_per,v_per
    # In[]
def sat_ECI(eccentricity,ecc_anomaly,a_semi_major_axis,omega_longitude_ascending_node,omega_argument_periapsis,inclination,nt_mean_motion):
    #Perifocal to ECI
    #Earth Centred Inertial Frame
    
    #Finds Perifocal Components
    r_per,v_per=perifocal(eccentricity,ecc_anomaly,a_semi_major_axis,omega_longitude_ascending_node,omega_argument_periapsis,inclination,nt_mean_motion)
    
    
    #Creates transformation
    Per_to_ECI=R.from_euler('ZXZ',[-float(omega_longitude_ascending_node),-float(inclination),-float(omega_argument_periapsis)],degrees=True)
    #Note: RAAN,omega,inc in degrees
    pos_ECI=(Per_to_ECI.apply(r_per)).tolist()
    vel_ECI=(Per_to_ECI.apply(v_per)).tolist()
    
    
    
    
    
    return pos_ECI,vel_ECI
    # In[]
def sat_ECF(theta_t,eci_position,eci_velocity):
    
    
    #Creates rotational transformation
    ECI_to_ECF=R.from_euler('Z',[float(theta_t)], degrees=True)
    
    #Applies Rotational Transformation
    pos_ECF=ECI_to_ECF.apply(eci_position).tolist()[0] # had to perform weird 
    #weird conversion to get back to list
    
    vel_ECF=ECI_to_ECF.apply(eci_velocity).tolist()[0]
    #km/s
    
    
    #Relative Velocity
    Siderial_rotation=[1,1,360/86164.091] #Degrees/s
    vel_rel=ECI_to_ECF.apply(eci_velocity-np.matmul(Siderial_rotation,eci_position))
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
    
    station_longitude=float(station_longitude)*math.pi/180
    station_latitude=float(station_latitude)*math.pi/180
    #input as Rads only
    
    
    
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
    
    v_rel=range_topo_velocity
    
    #Calculates the AZ and EL
    AZ=math.atan(R[0]/R[1])
    EL=math.atan(R[2]/(math.sqrt(R[0]**2+R[1]**2)))
    
    
    r=np.linalg.norm(R) #scalar of R
    
    R_xy=[R[0],R[1]]
    #In Software Specification Rxy is [tx ty]{Rtx;Rty} which would give a result of a singular value
    #Here we assume the Professor meant R_xy= the x and y components of R
    
    v_xy=[v_rel[0],v_rel[1]]
    
    #Calculates rates of AZ and EL
    rate_of_AZ=np.cross(v_xy,R_xy)
    rate_of_EL=(r*v_rel[2]-R[2]*np.dot(R_xy,v_xy)/r)/(r**2)
    
    return AZ,EL,rate_of_AZ,rate_of_EL

    # In[]


 ##                         Initialization Code
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
        self.argper=line2[34:42]
        self.meanan=line2[43:51]
        self.meanmo=line2[52:63]
        self.ndot=line1[33:43]
        self.n2dot=line1[44:50]
        self.bstar=line1[53:61]
        self.orbitnum=line2[63:68]
    # In[]
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
            print('Satellite has been registered as Satellite:',len(SatList),"Array index: [",len(SatList)-1,"]")
        except:
            SatList=[]
            SatList.append(Satellite(line0,line1,line2))
            print('SatList has been created')
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
def User_Input_parser_Call(StationLocationStr,TLEfilestr,TrackingSchedulestr,LinkInputsstr):
    StationInstance=Station(StationLocationStr)
    SatList=SatListPropagate(TLEfilestr)
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
        #degrees,
        ecc_anomaly=KeplerEqn(Mt_Mean_anomaly,SatList[p].eccn)
        
        
        #returns in radians
        mu=398600.4418 #km^3/s^2
        a=(mu/(2*np.pi*float(SatList[p].meanmo)/86400)**2)**(1/3)
        [pos_ECI,vel_ECI]=sat_ECI(SatList[p].eccn,KeplerEqn(SatList[p].meanan,SatList[p].eccn), \
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
        
        zTest_ECI_R.append(pos_ECI),zTest_ECI_R.append(vel_ECI)
        zTest_ECF_R.append(pos_ECF),zTest_ECF_vel.append(vel_ECF),zTest_ECF_vel_rel.append(vel_rel_ECF)
        zTest_T_x.append(Tx),zTest_T_y.append(Ty),zTest_T_z.append(Tz)
        
        
        
        
        
      Time_dt=Time_dt+dt.timedelta(seconds=float(Tracking.timestep))
      #At the End change Time
      
    
    
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
  
  

  
  for i in range(0,int(StnInstance.az_el_nlim)):
      
      # For Each iteration of the Station Instation Limits
    ThisStationLimit=StnInstance.az_el_lim[i].split(",")
    
    # Covnerts Limits from Degrees to Rads
    AZ_muth_limit=float(ThisStationLimit[0])*math.pi/180
    EL_lim_max=float(ThisStationLimit[2])*math.pi/180
    EL_lim_min=float(ThisStationLimit[1])*math.pi/180
    
    
    for j in range(0,len(AZ_list)):
      
      if AZ_list[j] > (AZ_muth_limit) and (EL_lim_max) > EL_list[j] and EL_list[j] > (EL_lim_min):
        #compares Azimuth and Elevations to Limits
          
          #Satellite Available
          
        
          #Compares to Previous Value
          # If previous value is not in limits but is now either AOS or from another limit iteration** 
        if j >= (Satnum_iteration) and  Satnum_list[j] == Satnum_list[j-Satnum_iteration] and AZ_list[j-Satnum_iteration] < float(AZ_muth_limit) or float(EL_lim_max) < EL_list[j-Satnum_iteration] or EL_list[j-Satnum_iteration] < float(EL_lim_min):
            # A Signal has been acquired
            
            #** Here we test to see if the AOS was caused by limit bound issues or from actual AOS
            
            if i > 0:
                
                for k in range(0,len(Times_LOS)):
                    
                    if time[j] == Times_LOS[k] and Satnum[j]== SatNum_LOS[k]:
                        append=0
                        del Times_LOS[k]
                        del SatNum_LOS[k]
                        del AZ_LOS[k]
                        del EL_LOS[k]
                    else: 
                        append=1
                        break
                if append==1:
                    AZ_AOS.append(AZ_list[j])
                    EL_AOS.append(EL_list[j])
                    Times_AOS.append(time[j])
                    SatNum_AOS.append(Satnum_list[j])
                    Signal_lost_AOS.append(Signal_lost[j])
                        
                #Then no AOS, LOS
                
            # This is actual AOS    
            else:
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
      else:
          # Satellite Unavailable
          
          #Compares previous iteration of Satellite to current Satellite Available
          if j >= (Satnum_iteration)  and Satnum_list[j]==Satnum_list[j-(Satnum_iteration)] and AZ_list[j-(Satnum_iteration)] > float(AZ_muth_limit) and float(EL_lim_max) > EL_list[j-(Satnum_iteration)] and EL_list[j-(Satnum_iteration)] > float(EL_lim_min):
              #If Satellite was available but is no longer then it is either LOS or out of our Limit Range
              #
              
              AZ_LOS.append(AZ_list[j])
              EL_LOS.append(EL_list[j])
              Times_LOS.append(time[j])
              SatNum_LOS.append(Satnum_list[j])
              Signal_lost_LOS.append(Signal_lost[j])
              
              

    
  AOS_List=[AZ_AOS,EL_AOS,Times_AOS,SatNum_AOS,Signal_lost_AOS]
  LOS_List=[AZ_LOS,EL_LOS,Times_LOS,SatNum_LOS,Signal_lost_LOS]
#Creates a list of Available Azimuth, Elevation, Times and Satnum available 
  
  
  return AZ_avail,EL_avail,Times_avail,Satnum_avail,AOS_List,LOS_List


    # In[]
#def Link_Calculations(LinkData):


def linkcal(linkdat):
  Linkcalcfile=open(linkdat, 'rt')
  frequency=Linkcalcfile.readline()
  Antennaeff=Linkcalcfile.readline()
  AntennaDia=Linkcalcfile.readline()
  Linkcalcfile.close
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
            if Satnum_AOS[i] == Satnum_LOS[j]:
                Templist=[Satnum_AOS[i],SatList[Satnum_AOS[i]].name,Sat_AOS_Time[i],Sat_LOS_Time[j],Sat_Signal_Lost[i]]
                #creates a temporary list to store Sat number,Sat list name ,AOS time,LOS time and 
                AOS_LOS_list.append(Templist)
                
        
    
    
    
    
    return AOS_LOS_list
    # In[]
def TrackingData(freq,Antennaeff,AntennaDia,R_ti):
    Signal_loss=[]
    for i in range(0,len(R_ti)):
        R=np.linalg.norm(R_ti[i])
        Signal_loss.append(20*math.log(10,4*math.pi*(float(R))/(3.0e8*float(freq)*1e6)))
    return Signal_loss

 
   

    # In[]
#       Debugging 
    
    
    # In[]
def STKout(EphemFile,StartString,time,Coord,position,velocity):
    
    #Asssuming StartString is the Start Time String
    
    #Here I interpreted the positon and velocity arrays as
    # Postion Array:
    #  X Y Z
    
    #opens or creates a file with the name EphemFile
    file_1=open(EphemFile,"w+")
    
    #Creates a set of Strings that can be easily changed. 
    
    ScenarioEpochString="ScenarioEpoch \t"+StartString
    CentralBody="Earth"
    CentralBodyString="Central Body "+ CentralBody
    CoordinateSystemString="CoordinateSystem "+Coord
    NumofEphermis=len(time) #I don't know how to attain this number
    
    #Writing Header
    file_1.write("stk.v.4.3 \n\nBEGIN Ephemeris \n\n")
    file_1.write("NumberOfEpemerisPoints "+str(NumofEphermis)+"\n")
    file_1.write(ScenarioEpochString+"\n")
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
    while i<len(time):
        
        file_1.write('{0:15.14e} {1:15.14e} {2:15.14e} {3:15.14e} {4:15.14e} {6:15.14e} \n'.format(float(time[i]),float(X[i]),float(Y[i]),float(Z[i]),float(X_dot[i]),float(Y_dot[i]),float(Z_dot[i])))
        #writing to file in formatted way
        
        
        
        # for i is in time, we iterate through the list to write in the values to the value
        i=i+1
        
    file_1.write("\n\nEND Ephemeris")
    file_1.close()

# In[]
def STKsp(PointingAngles,time):
        
        #We should have an input file name. I've provided one here to help
        spFile="STKsp.sp"
        
        #opens file
        file_2=open(spFile,"w+")
        
        NumofAttitudes=len(time) 
        
        #Writing Header
        file_2.write("stk.v.4.3 \nBEGIN\tAttitude\n")
        file_2.write("NumberOfAttitudePoints "+str(NumofAttitudes)+"\n")
        file_2.write("AttitudeTimeAzElAngles\n")
        AL=[row[0] for row in PointingAngles]
        EL=[row[1] for row in PointingAngles]
        i=0
        while i<len(time):
            file_2.write('{0:7.2f} {1:7.2f} {2:7.2f}\n'.format(float(time[i]),float(AL[i]),float(EL[i])))
        # for i is in time, we iterate through the list to write in the values to the value
            i=i+1
            
        #writing closer 
        file_2.write("END Attitude")
        file_2.close    

 # In[]

##                  Main Function


#Assuming That The Use has alreadu initialized all necessary functions and classes
# The Main Program can be deduced to this
[StationInstance,SatList,Tracking,LinkData]=User_Input_parser_Call(r'D:\School\5th Year Fall Semester\ESSE 4350\Lab 03\ReferenceFiles\Station.txt',r'D:\School\5th Year Fall Semester\ESSE 4350\MiniQuiz2\gps-ops.txt',r'D:\School\5th Year Fall Semester\ESSE 4350\MiniQuiz2\TrackingData.txt',r'D:\School\5th Year Fall Semester\ESSE 4350\Lab 03\ReferenceFiles\LinkInputs.txt')
[AZ,EL,Rate_of_AZ,Rate_of_EL,R_ti,v_rel_ti,time,Satnum]=Sat_pos_velCall(StationInstance,SatList,Tracking)

[freq,Antennaeff,AntennaDia]=linkcal(r'D:\School\5th Year Fall Semester\ESSE 4350\Lab 03\ReferenceFiles\LinkInputs.txt')
Signal_loss=TrackingData(freq,Antennaeff,AntennaDia,R_ti)
[AZ_avail,EL_avail,Times_avail,Satnum_avail,AOS_List,LOS_List]=Pointing(StationInstance,AZ,EL,time,Satnum,Signal_loss)
#Visibility creates a formatted list 
[AOS_LOS_list]=Visibility(StationInstance,AZ,EL,time,Satnum,Signal_loss)

# Debugging
STKout('EphemFileExample.txt',str(Starttime),time,'Inertial',zTest_ECI_R,zTest_ECI_v)



#Outputs AZ in Rads

    # In[]
signalloss=linkcal(r'D:\School\5th Year Fall Semester\ESSE 4350\Lab 03\ReferenceFiles\LinkInputs.txt')

#Note: AOS/LOS has not been created yet

    # In[]
#                                           Testing Cell

Refepoch=referenceepoch_propagate(Tracking)
#THETAN has been edited to input time_start_dt and Time_dt
#Tracking is Trackinginstance data
Starttime=dt.datetime.strptime(Tracking.starttime,'%Y-%m-%d-%H:%M:%S')
Midtime=Starttime.replace(hour=0,minute=0,second=0)
J2000=dt.datetime.strptime('2000-01-01 12:00:00','%Y-%m-%d %H:%M:%S')
D_u=(Midtime-J2000).days+(Midtime-J2000).seconds/86400
T_u=D_u/36252
D_u_2=(Starttime-J2000).days+(Starttime-J2000).seconds/86400
T_u_2=D_u_2/36525
GMST_00=(99.9677947+36000.7700631*T_u+0.00038793*T_u**2-2.6e-8*T_u**3)%360
r=1.002737909350795+5.9006e-11*(T_u)-5.9e-15*(T_u)**2
del_sec=(Starttime-Midtime).total_seconds()
GMST_t=(GMST_00+360*r/86400*(del_sec))%360

GMST=THETAN(Refepoch)
Time_dt=dt.datetime.strptime(Tracking.starttime,'%Y-%m-%d-%H:%M:%S')
p=17
[Mt_Mean_anomaly,Nt_anomaly_motion]=mean_anomaly_motion(Time_dt,SatList[p].refepoch,float(SatList[p].meanan),float(SatList[p].meanmo),float(SatList[p].ndot),float(SatList[p].n2dot))
        #degrees,
ecc_anomaly=KeplerEqn(Mt_Mean_anomaly,SatList[p].eccn)
        #returns in radians
mu=398600.4418 #km^3/s^2
a=(mu/(2*np.pi*float(SatList[p].meanmo)/86400)**2)**(1/3)
[pos_ECI,vel_ECI]=sat_ECI(SatList[p].eccn,KeplerEqn(SatList[p].meanan,SatList[p].eccn), \
        a,SatList[p].raan,SatList[p].argper,SatList[p].incl,Nt_anomaly_motion)