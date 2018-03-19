#!/bin/bash

. ../SETTINGS.sh

## PART III
## ==============================================

# now push release branch
cd SGpp
git push

# push Github release repository
cd ../SGpp_github/
git push
cd ..
