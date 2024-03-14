mkdir test_results
"%WINSGPP_PATH%\bin\base\Debug\tests.exe" --log_level=unit_scope --log_sink=test_results\base.log
"%WINSGPP_PATH%\bin\solver\Debug\tests.exe" --log_level=unit_scope --log_sink=test_results\solver.log
"%WINSGPP_PATH%\bin\combigrid\Debug\tests.exe" --log_level=unit_scope --log_sink=test_results\combigrid.log
"%WINSGPP_PATH%\bin\datadriven\Debug\tests.exe" --log_level=unit_scope --log_sink=test_results\datadriven.log --run_test=!DBMatDatabaseTest/TestOverwriteMatrix:!HPOTest/fitScalesGP:!HPOTest/harmonicaConfigs
"%WINSGPP_PATH%\bin\optimization\Debug\tests.exe" --log_level=unit_scope --log_sink=test_results\optimization.log
"%WINSGPP_PATH%\bin\pde\Debug\tests.exe" --log_level=unit_scope --log_sink=test_results\pde.log
"%WINSGPP_PATH%\bin\quadrature\Debug\tests.exe" --log_level=unit_scope --log_sink=test_results\quadrature.log