This is a recipe for a conan package. Simply run

conan create . sgpp/experiment

to create the package! It downloads an builds SGpp (3.2.0) in the conan directory. It is then available to use as a normal conan package.
For example, a application using SGpp, fmt and ranges-v3 would then be able to express the dependencies like this:

[requires]
fmt/6.0.0@bincrafters/stable
range-v3/0.10.0@ericniebler/stable
SGpp/3.2.0@sgpp/experimental

[generators]
cmake

This would generate a cmake fragment (once the user calls "conan install <file>") that handles all dependencies listed above! Note that it is also possible to create a scons fragment (see generator).
