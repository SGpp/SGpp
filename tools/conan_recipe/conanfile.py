from conans import ConanFile, CMake, tools


class SgppConan(ConanFile):
    name = "SGpp"
    version = "3.2.0"
    license = "BSD-style license"
    author = "Gregor Daiß Gregor.Daiss@ipvs.uni-stuttgart.de"
    url = "https://sgpp.sparsegrids.org/"
    description = """The sparse grids toolkit SG++
 SG++ is a collection of numerical algorithms for sparse grids. It
 contains modules for interpolation, quadrature, data mining
 (regression, classification, clustering), optimization, PDEs, and
 more. SG++ implements algorithms for spatially adaptive grids and
 also provides a module for the combination technique. Many of the
 implemented algorithms are also available as a high-performance
 version, often orders of magnitude faster than standard
 implementations."""
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = {"shared": False}
    generators = "cmake"

    def source(self):
        self.run("git clone https://github.com/SGpp/SGpp")
        with tools.chdir("{}/SGpp".format(self.source_folder)):
            self.run("pwd")
            self.run("git checkout v3.2.0")

    def build(self):
        with tools.chdir("{}/SGpp".format(self.source_folder)):
            self.run('scons -j4 SG_JAVA=0 RUN_BOOST_TESTS=0')

    def package(self):
        self.copy("*.hpp", dst="include")
        self.copy("*.so", dst="lib", keep_path=False)

    def package_info(self):
        #self.cpp_info.libs = ["sgppbase"]
        self.cpp_info.libs = ["sgppbase", "sgppquadrature", "sgppmisc", "sgppdatadriven", "sgpppde", "sgppsolver", "sgppoptimization", "sgppcombigrid"]

