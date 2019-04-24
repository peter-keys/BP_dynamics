# BP_dynamics
Some codes for looking at MBP dynamics (specifically waves)

I have put up some codes here for looking at MBP dynamics in data and si mulations. 
These are specifically useful for looking at wave motions in MBPs.

The codes are separtated into different sub-repositories depending on their use. 
The Sub-repositories are:

Inversion_work: 
These are codes specifically designed to return useful information from inverison codes on MBPs.
That is, these codes look at preparing MBP data for inversion in NICOLE and SIR, as well as 
returning key values from the code such B-field etc, for these MBPs.
Also, there are some codes here that look at patching a temporal sub-field of the MBP as 
a 2D array to invert to make it easier to invert a sub field for a specific MBP.

Velocity_code:
This needs to be run in conjunction with the 'tracking_code' repository. This finds the 
velocity characteristics of a tracked field of view. This can be used with codes in the 
'Inversion_work' subdirectory to get the inverted characteristics for the MBPs in the same file 
for temporal comparsions. This is used to establish MBPs that experience motions greater than 2km/s 
for driving kink oscillations.

Vortex_work:
This is a work in progress. This will be a code that finds vortex flows within the field-of-view. 
This can then be compared to the tracked MBPs to find those that fall into regions of higher vorticity 
This will then be used as a basis of detecting MBPs with possible Alfvenic motions. The basis of this 
work is too use Local Correlation Tracking to get the horizontal velocity component which will then be 
combined with teh line-of-sight velocity component to work out the 3D velocity component and therefore 
the vorticity. The code needs to have the LOS component added and needs to be tested against simulations 
with know vorticity. There are other codes within this repository which can be used to establish 
supergranule boundaries from HMI continuum images.

tracking_code:
This detects and tracks MBP locations in an intensity image array. This is then used to find candidates 
for waves in MBPs. 

Wave_code:
