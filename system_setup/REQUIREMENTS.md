# Requirements setup for HistoneDB 

To fully and correctly run the project, some system configurations and dependencies are important. To setup Ubuntu 18.04 machine system, python packages and MySQL server run script ``` system_setup/full_setup.sh ```:
```
sh system_setup/full_setup.sh
```
This command do: 

1) Install all system dependencies (see ```system_setup/sys_requirements_setup.txt```) and configure MySQL server

2) Create and activate virtualenv:

```
virtualenv histdb_py27
source histdb_py27/bin/activate
```

3) Install all python dependencies (see ```system_setup/py_requirements_setup.txt```).

***Note***: the step includes installing mysql-client (without password).

