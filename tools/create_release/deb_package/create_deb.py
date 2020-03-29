#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# TODO These packages do not really support changelogs yet - This should be added!

import sys
import os
import re
import shutil
import subprocess
import jinja2
#Environment, PackageLoader, select_autoescape
jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader("templates"))

import copy_headers

print("Welcome to the SG++ deb-package creation utility")
print("(please run this utility from the \"deb_pkg\" folder)")

# make sure that SG++ was actually built
if not os.path.exists("../../../lib/libsgppbase.so"):
    print("Aborting: At least \"libsgppbase.so\" has to exist in \"../../../lib\"")
    sys.exit(1)

# make sure that SG++ was built with python spport
if not os.path.exists("../../../lib/pysgpp/_pysgpp_swig.so"):
    print("Aborting: \"_pysgpp_swig.so\" not found, SG++ was probably not build with python support")
    sys.exit(1)

# figure out which modules were built
modules = []
for so_file in os.listdir("../../../lib/"):
    if not so_file.endswith(".so"):
        continue
    modules.append(re.search(r"libsgpp(.*?)\.so", so_file).group(1))
print(("Found compiled modules: " + str(modules)))

# check for dependencies
print("Building dependencies")
sgpp_deps = ""
pysgpp_deps = "libsgpp"
jsgpp_deps = "libsgpp"
if len(sys.argv) >= 2:
    with open(sys.argv[1], 'r') as fin:
        line = fin.readline()
        while line:
            sgpp_deps += (", " + line.strip())
            print("\t" + line.strip())
            line = fin.readline()
    sgpp_deps = sgpp_deps.strip(" ,")
    if len(sys.argv) >= 3:
        with open(sys.argv[2], 'r') as fin:
            line = fin.readline()
            while line:
                pysgpp_deps += (", " + line.strip())
                print("\t" + line.strip())
                line = fin.readline()
        pysgpp_deps = pysgpp_deps.strip(" ,")
        if len(sys.argv) >= 4:
            with open(sys.argv[2], 'r') as fin:
                line = fin.readline()
                while line:
                    jsgpp_deps += (", " + line.strip())
                    print("\t" + line.strip())
                    line = fin.readline()
            jsgpp_deps = jsgpp_deps.strip(" ,")

package_name = ""
if len(sys.argv) >= 5:
    print("WARNING! Running with default version numbers for automated testing!")
    major_version = "0"
    minor_version = "0"
    package_revision = "0"
    maintainer_name = "dummy"
    maintainer_email = "dummy"
    package_name = "libsgpp-test-package_" + major_version + "." + minor_version + "-" + package_revision
else:
    # enter version information
    major_version = input("Enter major version: ")
    if major_version == "":
        print("Aborting: You have to specify a major version")
        sys.exit(1)
    if not major_version.isdigit():
        print("Aborting: major version can only contain digits")
        sys.exit(1)
    minor_version = input("Enter minor version (default: \"0\"): ")
    if minor_version == "":
        minor_version = "0"
    if not minor_version.isdigit():
        print("Aborting: minor version can only contain digits")
        sys.exit(1)
    package_revision = input("Enter package_revision (default: \"1\"): ")
    if package_revision == "":
        package_revision = "1"
    if not package_revision.isdigit():
        print("Aborting: package revision can only contain digits")
        sys.exit(1)

    # enter maintainer information
    maintainer_name = input("Enter the name of the maintainer (default \"Dirk Pflüger\"):")
    if maintainer_name == "":
        maintainer_name = "Dirk Pflüger"

    maintainer_email = input("Enter the email adress of the maintainer (default \"Dirk.Pflueger@ipvs.uni-stuttgart.de\"):")
    if maintainer_email == "":
        maintainer_email = "Dirk.Pflueger@ipvs.uni-stuttgart.de"

    package_name = "libsgpp_" + major_version + "." + minor_version + "-" + package_revision
print(("package name: " + package_name))
deb_name = package_name + ".deb"

if os.path.exists(package_name):
    print(("Aborting: build path \"" + package_name + "\" already exists (remove or create different version)"))
    sys.exit(1)

if os.path.exists(deb_name):
    print(("Aborting: The deb file to build (\""+ deb_name + "\") already exists"))
    sys.exit(1)

# create basic package structure
os.makedirs(package_name + "/usr/lib")
os.makedirs(package_name + "/usr/include")
# copy library files
for module in modules:
    shutil.copy("../../../lib/libsgpp" + module +".so", package_name + "/usr/lib")
# copy C++-headers
for module in modules:
    copy_headers.copy_headers(module, package_name)
# copy meta-files
os.makedirs(package_name + "/DEBIAN")
control_template = jinja2_env.get_template('control_template')
rendered = control_template.render(major_version=major_version, minor_version=minor_version, package_revision=package_revision, maintainer_name=maintainer_name, maintainer_email=maintainer_email, sgpp_deps=sgpp_deps)
with open(package_name + "/DEBIAN/control", "w") as f:
    # dpkg-deb requires newline at the end of the file
    f.write(rendered + "\n")

# finally, create the package
try:
    subprocess.check_call(["fakeroot", "dpkg-deb", "--build", package_name])
except Exception as e:
    print("Error: building the deb-file failed")
    print("reason: ")
    print((str(e)))

# now create the python bindings package
if len(sys.argv) >= 5:
    package_name = "libsgpp-python-test-package_" + major_version + "." + minor_version + "-" + package_revision
else:
    package_name = "libsgpp-python_" + major_version + "." + minor_version + "-" + package_revision
print(("python bindings package name: " + package_name))
deb_name = package_name + ".deb"

if os.path.exists(package_name):
    print(("Aborting: build path \"" + package_name + "\" already exists (remove or create different version)"))
    sys.exit(1)

if os.path.exists(deb_name):
    print(("Aborting: The deb file to build (\""+ deb_name + "\") already exists"))
    sys.exit(1)

os.makedirs(os.path.join(package_name, "usr/lib/python3/dist-packages/pysgpp/"))

# copy main shared library
shutil.copy("../../../lib/pysgpp/_pysgpp_swig.so", os.path.join(package_name, "usr/lib/python3/dist-packages/pysgpp/_pysgpp_swig.so"))
# copy library files
import copy_python
copy_python.copy_python(package_name)
# copy meta-files
os.makedirs(os.path.join(package_name, "DEBIAN"))
control_template = jinja2_env.get_template('control_python_template')
rendered = control_template.render(major_version=major_version, minor_version=minor_version, package_revision=package_revision, maintainer_name=maintainer_name, maintainer_email=maintainer_email, pysgpp_deps=pysgpp_deps)
with open(os.path.join(package_name, "DEBIAN/control"), "w") as f:
    # dpkg-deb requires newline at the end of the file
    f.write(rendered + "\n")

# finally, create the package
try:
    subprocess.check_call(["fakeroot", "dpkg-deb", "--build", package_name])
except Exception as e:
    print("Error: building the deb-file failed")
    print("reason: ")
    print((str(e)))

# now create the java bindings package
if len(sys.argv) >= 5:
    package_name = "libsgpp-java-test-package_" + major_version + "." + minor_version + "-" + package_revision
else:
    package_name = "libsgpp-java_" + major_version + "." + minor_version + "-" + package_revision
print(("java bindings package name: " + package_name))
deb_name = package_name + ".deb"

if os.path.exists(package_name):
    print(("Aborting: build path \"" + package_name + "\" already exists (remove or create different version)"))
    sys.exit(1)

if os.path.exists(deb_name):
    print(("Aborting: The deb file to build (\""+ deb_name + "\") already exists"))
    sys.exit(1)

os.makedirs(os.path.join(package_name, "debian/usr/share/java/"))
os.makedirs(os.path.join(package_name, "debian/usr/lib/"))
#os.makedirs(os.path.join(package_name, "DEBIAN"))

# copy main shared library
shutil.copy("../../../lib/jsgpp/libjsgpp.so", os.path.join(package_name, "debian/usr/lib/libjsgpp.so"))

control_template = jinja2_env.get_template('control_java_template')
rendered = control_template.render(major_version=major_version, minor_version=minor_version, package_revision=package_revision, maintainer_name=maintainer_name, maintainer_email=maintainer_email, jsgpp_deps=jsgpp_deps)
with open(os.path.join(package_name, "debian/control"), "w") as f:
    # dpkg-deb quires newline at the end of the file
    f.write(rendered + "\n")

changelog_template = jinja2_env.get_template('changelog_java_template')
rendered = changelog_template.render(major_version=major_version, minor_version=minor_version, package_revision=package_revision)
with open(os.path.join(package_name, "debian/changelog"), "w") as f:
    # dpkg-deb quires newline at the end of the file
    f.write(rendered + "\n")
# Run javahelper! For some reason java helper expects the DEBIAN folder to be names debian and installs everything in there?
# This is rather odd: Maybe it's not meant to be used like this - it works however. 
# TODO Find alternative to javahelper (or at least use it in a way to get rid of the following moves)
os.chdir(package_name)
os.system("jh_installlibs ../../../../lib/jsgpp/jsgpp.jar")
os.system("mv debian/usr usr")
os.system("mv debian DEBIAN")
os.chdir("..")
# finally, create the package
try:
    subprocess.check_call(["fakeroot", "dpkg-deb", "--build", package_name])
except Exception as e:
    print("Error: building the deb-file failed")
    print("reason: ")
    print((str(e)))
