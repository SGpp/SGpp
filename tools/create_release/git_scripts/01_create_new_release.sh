#!/bin/bash

. ../SETTINGS.sh

## PART I
## ==============================================

# clone Github release repository
git clone git@github.com:SGpp/SGpp.git SGpp_github

# clone SG++ without checkout
git clone --no-checkout git@simsgs.informatik.uni-stuttgart.de:SGpp/SGpp.git
# checkout based on whitelist
cd SGpp
(git checkout ${REVISION} -- `cat ../../RELEASE_CONTENTS.sh | grep -v "#"`) || exit 1
# create new release branch
git checkout -b ${RELEASE}
git commit -am "New release branch ${RELEASE}"

# write new state for Github release repository
rsync --exclude=.git -av --delete-after . ../SGpp_github/
cd ../SGpp_github/
git add .
# checkin
git commit -m "New release ${RELEASE}"
cd ..

echo "Now test!"
