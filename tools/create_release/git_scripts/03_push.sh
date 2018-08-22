#!/bin/bash

. ../SETTINGS.sh

## PART III
## ==============================================

# now push release branch
cd SGpp
git push --set-upstream origin ${RELEASE}

# push Github release repository
cd ../SGpp_github/
git push
git push origin ${RELEASE}
cd ..
