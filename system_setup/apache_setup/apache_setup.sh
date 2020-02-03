#!/bin/bash

parentdir="$(dirname "$PWD")"
parentdir="$(dirname "$parentdir")"
sed -i 's#MY_PROJECT_DIRECTORY#'$parentdir'#g' wsgi.conf

sudo cp -f wsgi.conf /etc/apache2/mods-available