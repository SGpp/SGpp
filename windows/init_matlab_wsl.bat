@echo off

if [%1]==[] goto usage
wsl sudo ./matlab_wsl/requirements.sh;sudo ./matlab_wsl/install.sh %1;
goto :eof
:usage
@echo Usage: %0 ^<matlab download^>
pause