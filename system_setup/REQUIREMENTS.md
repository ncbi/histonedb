# Requirements setup for HistoneDB 

To fully and correctly run the project, some system configurations and dependencies are important. To setup Ubuntu machine system, python packages and MySQL database follow the steps below:

1) Create a virtualenv (for correct setup, please, do not use anaconda) in your working directory and activate it
```
virtualenv histdb_py27
source histdb_py27/bin/activate
``` 

2) For system setup run script system_setup/full_setup.sh
```
sh system_setup/full_setup.sh
```
This command will install all system dependencies (see ```system_setup/sys_requirements_setup.txt```) and all python dependencies in activated virtualenv (see ```system_setup/py_requirements_setup.txt```). ***Note***: the step includes installing mysql-client.

