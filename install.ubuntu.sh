#!/bin/bash
#
# This script is designed install pre-requisites for PrimerSeq on linux
# (Debian based distros with apt-get). The following is assumed when you run
# this script:
#
# You have root access
# The default python version is 2.7.X

EXPECTED_ARGS=0
EXPECTED_ARGS=65

if [ $# -ne $EXPECTED_ARGS ]
then
	echo "This script is designed to install pre-requisites for PrimerSeq on linux (Debian based distros with apt-get). It is assumed you have the following:\n* super user access\n* Default Python version is 2.7"
	exit $E_BADARGS
fi

# Get dependencies of pip
sudo apt-get python-setuptools
sudo apt-get python-pip

# possible dependencies for matplotlib
sudo apt-get install libtiff
sudo apt-get install libpng

# install python dependencies
sudo pip install numpy
sudo pip install matplotlib
sudo pip install pygr
sudo pip install networkx

# add wxwidgets to sources list
deb http://apt.wxwidgets.org/ `lsb_releas -a | grep Codename | cut -f2`-wx main
deb-src http://apt.wxwidgets.org/ `lsb_releas -a | grep Codename | cut -f2`-wx main

# update so wxwidgets repo is found with apt-get
sudo apt-get update

# install wxPython
sudo apt-get install python-wxgtk2.8 python-wxtools wx2.8-il8n libwxgtk2.8-dev libgtk2.0-dev