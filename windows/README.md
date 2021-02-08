## Compiling and debugging sgpp on Windows using Visual Studio

- Install [Visual Studio 19](https://visualstudio.microsoft.com/). This is totally different from Vs Code!
- Set the environment variable `WINSGPP_PATH` to the cloned repository directory
- Add `%WINSGPP_PATH%\lib\Debug` to your `PATH` variable if you want to be able to run examples or tests
- Open the solution file `winsgpp.sln` located in the root directory with Visual Studio
- You can then build and run SGpp in Visual Studio
- Alternatively, you can execute the scripts `build_msbuild.bat` to build and `run_tests.bat` to run the tests.
  Note that you need to install the [VS build tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019) for that
