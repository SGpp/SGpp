#!/bin/bash

mkdir -p collectedStatistics
scp pfandedd@vgpu2:/scratch/pfandedd/git/SGppAutoTune/statistics/* collectedStatistics/
scp pfandedd@neon:~/git/SGppAutoTuneNeon/statistics/* collectedStatistics/
scp pfandedd@kepler:~/git/SGppAutoTuneKepler/statistics/* collectedStatistics/
scp pfandedd@pcsgs07:~/git/SGppAutoTune/statistics/* collectedStatistics/
scp winter@itojvkowjmf4u8jk.myfritz.net:~/SGpp/statistics/* collectedStatistics/

mkdir -p collectedLogs
scp pfandedd@vgpu2:/scratch/pfandedd/git/SGppAutoTune/vgpu2.log collectedLogs/
scp pfandedd@neon:~/git/SGppAutoTuneNeon/neon.log collectedLogs/
scp pfandedd@kepler:~/git/SGppAutoTuneKepler/kepler.log collectedLogs/
scp pfandedd@pcsgs07:~/git/SGppAutoTune/tahiti.log collectedLogs/
scp winter@itojvkowjmf4u8jk.myfritz.net:~/SGpp/hawaii.log collectedLogs/
