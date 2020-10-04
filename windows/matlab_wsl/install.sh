export DISPLAY=":0"
unzip $(wslpath $1) -d /tmp/matlab_install/
cd /tmp/matlab_install/bin/glnxa64/
sudo ./install_unix_legacy