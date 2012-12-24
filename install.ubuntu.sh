#!/bin/sh
#
# This script is designed install pre-requisites for PrimerSeq on linux
# (Debian based distros with apt-get). The following is assumed when you run
# this script:
#
# You have root access
# The default python version is 2.7.X

tty -s; if [ $? -ne 0 ]; then konsole -e "su root - -c $0"; exit; fi
EXPECTED_ARGS=0
E_BADARGS=65

if [ $? -ne $EXPECTED_ARGS ]
then
	echo "This script is designed to install pre-requisites for PrimerSeq on linux (Debian based distros with apt-get). It is assumed you have the following:\n* root access\n* Default Python version is 2.7"
	exit $E_BADARGS
fi

# Get dependencies of pip
sudo apt-get install python-setuptools
sudo apt-get install python-pip

# possible dependencies for matplotlib
sudo apt-get install libtiff
sudo apt-get install libpng

# install python dependencies
sudo pip install numpy
sudo pip install matplotlib
sudo pip install pygr
sudo pip install networkx

# add wxwidgets to sources list
sudo echo "deb http://apt.wxwidgets.org/ `lsb_release -a | grep Codename | cut -f2`-wx main
deb-src http://apt.wxwidgets.org/ `lsb_release -a | grep Codename | cut -f2`-wx main" > /etc/apt/sources.list.d/wxwidgets.list

# update so wxwidgets repo is found with apt-get
sudo apt-get update

# install wxPython
sudo apt-get install python-wxgtk2.8 python-wxtools wx2.8-il8n libwxgtk2.8-dev libgtk2.0-dev

echo "You may want to copy the installation results. Press enter when finished . . ."
read -p "$done"
