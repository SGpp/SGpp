mkdir -p statistics
scp pfandedd@vgpu2:/scratch/pfandedd/git/SGppAutoTune/statistics/* statistics/
scp pfandedd@neon:~/git/SGppAutoTuneNeon/statistics/* statistics/
scp pfandedd@kepler:~/git/SGppAutoTuneKepler/statistics/* statistics/
scp pfandedd@pcsgs07:~/git/SGppAutoTune/statistics/* statistics/
scp winter@itojvkowjmf4u8jk.myfritz.net:~/git/SGppAutoTune/statistics/* statistics/

mkdir -p logs
scp pfandedd@vgpu2:/scratch/pfandedd/git/SGppAutoTune/vgpu2.log logs/
scp pfandedd@neon:~/git/SGppAutoTuneNeon/neon.log logs/
scp pfandedd@kepler:~/git/SGppAutoTuneKepler/kepler.log logs/
scp pfandedd@pcsgs07:~/git/SGppAutoTune/tahiti.log logs/
scp winter@itojvkowjmf4u8jk.myfritz.net:~/git/SGppAutoTune/hawaii.log logs/
