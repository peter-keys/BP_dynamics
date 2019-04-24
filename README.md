# BP_dynamics
Some codes for looking at MBP dynamics (specifically waves)

I have put up some codes here for looking at MBP dynamics in data and simulations. 

These are specifically useful for looking at wave motions in MBPs.

Please note, these codes may not all be easy to follow/understand. I am working on making a version of these codes that is more user-friendly and is much easier to use for the general user. The ultimate idea is to have some applet that you can supply your data to, which then will allow you to determine candidates for various waver modes in MBPs (note this may take a while, though the basis of such a code is here alreasy, it just needs to be tied together and tested). 

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
This is various useful wave codes by myself (and a few by some others) which is implemented in the wave studies. 
The codes have a variation on the wavelet procedure in IDL, a way to establish the spatial power from wavelets 
and the empirical mode decomposition technique which is often used in tandem with wavelets. Also, described 
here are methods for finding sausge modes in pores, and how to estanblish theoretically whether you have a 
surface or body mode. Also, there is a code that makes an estimate on energy values for known oscillations 
based on a theoretical description.
