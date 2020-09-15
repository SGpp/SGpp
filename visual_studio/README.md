## Compiling and debugging sgpp on Windows using Visual Studio

- Install [Visual Studio 19](https://visualstudio.microsoft.com/) and include the workload `Linux development with C++`. This is totally different from Vs Code!
- Install the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) (WSL 1, not WSL 2!). This README assumes that you use Ubuntu
- Install the following packages on the WSL: `sudo apt install g++ gdb make rsync zip`
- Start ssh on the WSL using `sudo service ssh start`
- Fetch the remote headers of the WSL as described [here](https://devblogs.microsoft.com/cppblog/intellisense-for-remote-linux-headers/)
- Remove the `.default` prefix from the files `sgpp.vcxproj.user.default`, `visual_studio\linux.props.default`, `visual_studio\build_release.py.default`, `visual_studio\build_debug.py.default`
- Edit the `visual_studio\linux.props` according to the comments
- Open the solution file `sgpp.sln` with Visual Studio
- Build the project and run. If any errors occur, please ensure that your `linux.props` is correct and otherwise report them
- There is no need to update your environment variables on your WSL, this is done automatically when debugging with Visual Studio

#### Issue: Adding, deleting and renaming files within VS does not work

Adding, deleting and renaming files does not work inside of Visual Studio, since VS does not like wildcard inclusions.
Solution: Add, delete and rename project files in the Windows Explorer, click Unload Project and then click Reload Project.
A fix is currently in development.

#### TODO: Make bindings also compileable and debuggeable on the WSL
