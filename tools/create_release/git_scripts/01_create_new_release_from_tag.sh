#!/bin/bash

. ../SETTINGS.sh

## PART I
## ==============================================

# clone Github release repository
git clone git@github.com:SGpp/SGpp.git SGpp_github

# clone SG++ for certain release tag
git clone -b ${RELEASE} git@simsgs.informatik.uni-stuttgart.de:SGpp/SGpp.git
cd SGpp
# # create new release branch
# git checkout -b ${RELEASE} ${REVISION}
# # checkout based on whitelist
# mv .git ../
# rm -rf * .*
# mv ../.git/ ./
# (git checkout ${RELEASE}) || exit 1
# git commit -am "New release branch ${RELEASE}"

# write new state for Github release repository
rsync --exclude=.git -av --delete-after . ../SGpp_github/
cd ../SGpp_github/
git add .
# checkin
git commit -m "New release ${RELEASE}"
git tag -a ${RELEASE}
cd ..

echo "Now test (at least SGpp)!"
