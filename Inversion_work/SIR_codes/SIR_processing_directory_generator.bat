;NAME: 
;	ROSA_directory_generator.bat
;  
; CATEGORY:
;   File I/O
;  	
;PURPOSE:
;	Creates the standard directory structures that we use with the ROSA pipeline
;	This is run from IDL command line to make things easier for the user
;	Assumes you are creating the directories in a UNIX/LINUX type system 
;	May not work, therefore, if you are trying to run it on a laptop
;	You can just manually create the necessary directory structures 
;
;  ** NOTE: YOU NEED TO EDIT SOME OF THIS FILE TO GET THE RIGHT DIRECTORIES/FILTERS **
;	
;
;REQUIREMENTS:
; 	Uses IDL programs: FILE_MKDIR (directory creator) & FILE_CHMOD (rights changer)
;	(P. Keys Feb 2019)
;-
;------------------------------------------------------------------------------------------------


directory_path = '/data/solarstore3/phk/data/27Jul2014/Inversions/SIR/AllMBPpix/'

FOR fr = 0, 88 DO BEGIN &$
	FILE_MKDIR,directory_path+'frame_'+STRING(fr,FORMAT='(I3.3)')+'/',/NOEXPAND_PATH &$
	FILE_MKDIR,directory_path+'frame_'+STRING(fr,FORMAT='(I3.3)')+'/observed_profiles/',/NOEXPAND_PATH &$
	FILE_MKDIR,directory_path+'frame_'+STRING(fr,FORMAT='(I3.3)')+'/models_output/',/NOEXPAND_PATH &$
	FILE_MKDIR,directory_path+'frame_'+STRING(fr,FORMAT='(I3.3)')+'/models_profiles/',/NOEXPAND_PATH &$
	FILE_MKDIR,directory_path+'frame_'+STRING(fr,FORMAT='(I3.3)')+'/models_errors/',/NOEXPAND_PATH &$
ENDFOR

; Later to move files from one place to another...
;FILE_MOVE, ['*.pro', 'makefile', 'mydata.dat'], 'BACKUP'

;FILE_MKDIR,raw_directory_path+'prespeckle/'+PI+Target1+Filter2,/NOEXPAND_PATH
;FILE_MKDIR,raw_directory_path+'prespeckle/'+PI+Target1+Filter3,/NOEXPAND_PATH
;FILE_MKDIR,raw_directory_path+'prespeckle/'+PI+Target1+Filter4,/NOEXPAND_PATH

; Moves files if needed to the right directory based on frame number...
directory_path = '/data/solarstore3/phk/data/27Jul2014/Inversions/SIR/AllMBPpix/'
FOR fr = 11, 88 DO BEGIN &$
	FILE_MOVE, ['Observed_profile*_fr'+STRING(fr,FORMAT='(I3.3)')+'*.per','Straylight*f'+STRING(fr,FORMAT='(I3.3)')+'*.per'], directory_path+'frame_'+STRING(fr,FORMAT='(I3.3)')+'/' &$
ENDFOR
