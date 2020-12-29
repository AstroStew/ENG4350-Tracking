# Welcome!
## Each branch of this Repository corresponds to a lab completed in for Space Hardware Class'(ENG 4350) Satellite Tracking Portion of the Class

The Code follows the Suggested Programming Structure (see Suggested_Programming_Structure.png)- this is also outlined in the Software Specifications document.

### Inputs
Inputs are in File called "Reference Files"
- Station Location
- TLE file
- Tracking Schedule
- Link Inputs (calcalated from Link Budget Equation)

### Outputs 
- Aqusition of Signal (AOS) List
- Loss of Signal (LOS) List
- Azimuth/Elevation (AZ/EL) Tables

### Known Issues:
#### Documentation Issues
- Some Files are not properly sorted
#### Code Issues
- Mean Anomally Off by 1 degree 
- Eccentricity Anomaly off by 1.5 degrees
- True Anomaly off by around 2 degrees
- Perifocal frame x-axis is off by 1000km and reversed signs when compared to STK
- Semi-major axis is off by 1 km

All Code Issues compared with STK plotted Values
Last Revised: 23:13 28-12-2020 
