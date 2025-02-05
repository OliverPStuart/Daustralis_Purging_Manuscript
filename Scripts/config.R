# R config file for capture paper scripts for thesis
# Contains variables for directory
# Two different home directories, one for personal laptop and one for work

# We check the system info to get which home directory use

env_info <- Sys.getenv(names=T)

if(env_info[attr(env_info,"names") == "LOGNAME"] == "cursorially"){
HOME_DIR="/home/cursorially/Daus_WGS_Paper"
} else {
HOME_DIR="/Volumes/Alter/Daus_WGS_Paper"
}

DATA_DIR=paste0(HOME_DIR,"/Data")
SCRIPT_DIR=paste0(HOME_DIR,"/Scripts")
FIGURE_DIR=paste0(HOME_DIR,"/Figures")
WORKING_DIR=paste0(HOME_DIR,"/Analysis")
REF_DIR=paste0(HOME_DIR,"/References")

# The relevant script will then define a working directory, usually the containing directory
