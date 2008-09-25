#!/bin/bash

#
# Setup cron job: crontab -e
#

cd ~/TUM/cronjob/SGpp

# check out / update
if [ -e ~/TUM/cronjob/SGpp/trunk ] 
then
    cd trunk
    svn up
else
    svn co file:///home_local/repositories/svn/SGpp/trunk
    cd trunk
fi

# compile
~/bin/doxygen

# copy
cp -r doc/ /home/wwwsccs/html/forschung/sg/SGpp/
chmod -R a+rx /home/wwwsccs/html/forschung/sg/SGpp/doc/
chmod -R g+w /home/wwwsccs/html/forschung/sg/SGpp/doc/

