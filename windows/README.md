## Compiling and debugging sgpp on Windows using Visual Studio

- Install [Visual Studio 19](https://visualstudio.microsoft.com/). This is totally different from Vs Code!
- Set the environment variable `WINSGPP_PATH` to the cloned repository directory
- Add `%WINSGPP_PATH%\lib\Debug` to your `Path` variable if you want to be able to run examples or tests
- Open the solution file `winsgpp.sln` located in the root directory with Visual Studio
- You can then build and run SGpp in Visual Studio

## Compiling sgpp on Windows using the MSBuild command line

Alternatively, you can also build and test SGpp just using the command line.
For that, you have to perform the following steps:

- Install the [VS build tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019)
- Define `WINSGPP_PATH` and add `%WINSGPP_PATH%\lib\Debug` to the `PATH` as described above
- Download [NuGet](https://www.nuget.org/downloads) and put `nuget.exe` either into the `PATH` or put it in the `SGpp/windows` directory
- Execute the scripts `build_msbuild.bat` to first build the library and then `run_tests.bat` to run the tests

## Compiling and running pysgpp on Windows

- Build the SGpp release configuration build as described in the previous section
- Define the environment variable `WINSGPP_PYTHON3` that points to the Python3 installation you want to use with pysgpp
- Build the project `pysgpp\pysgpp.vcxproj`

## TODO

- Change how source files are included. Use wildcards and wildcard exclusions
  of directories to make the build process easier when new source files are added
- Check all working directories of examples/tests, since all file read operations are using
  relative paths and executing them in the wrong working directory leads to errors
- Eliminate compiler warnings
- Fix failing tests

### Failing tests
- DBMatDatabaseTest/TestOverwriteMatrix
- HPOTest/fitScalesGP
- HPOTest/harmonicaConfigs
