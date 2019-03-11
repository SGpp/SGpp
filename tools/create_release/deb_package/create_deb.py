from builtins import input
#!/usr/bin/python3
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
if not os.path.exists("../../../lib/sgpp/libsgppbase.so"):
    print("Aborting: At least \"libsgppbase.so\" has to exist in \"../../../lib/sgpp\"")
    sys.exit(1)

# make sure that SG++ was built with python spport
if not os.path.exists("../../../lib/pysgpp/_pysgpp_swig.so"):
    print("Aborting: \"_pysgpp_swig.so\" not found, SG++ was probably not build with python support")
    sys.exit(1)

# figure out which modules were built
modules = []
for so_file in os.listdir("../../../lib/sgpp/"):
    if not so_file.endswith(".so"):
        continue
    modules.append(re.search(r"libsgpp(.*?)\.so", so_file).group(1))
print("Found compiled modules: " + str(modules))

# enter version information
major_version = str(eval(input("Enter major version: ")))
if major_version == "":
    print("Aborting: You have to specify a major version")
    sys.exit(1)
if not major_version.isdigit():
    print("Aborting: major version can only contain digits")
    sys.exit(1)
minor_version = str(eval(input("Enter minor version (default: \"0\"): ")))
if minor_version == "":
    minor_version = "0"
if not minor_version.isdigit():
    print("Aborting: minor version can only contain digits")
    sys.exit(1)
package_revision = eval(input("Enter package_revision (default: \"1\"): "))
if package_revision == "":
    package_revision = "1"
if not package_revision.isdigit():
    print("Aborting: package revision can only contain digits")
    sys.exit(1)

# enter maintainer information
maintainer_name = eval(input("Enter the name of the maintainer (default \"Dirk Pflüger\"):"))
if maintainer_name == "":
    maintainer_name = "Dirk Pflüger"

maintainer_email = eval(input("Enter the email adress of the maintainer (default \"Dirk.Pflueger@ipvs.uni-stuttgart.de\"):"))
if maintainer_email == "":
    maintainer_email = "Dirk.Pflueger@ipvs.uni-stuttgart.de"

package_name = "libsgpp_" + major_version + "." + minor_version + "-" + package_revision
print("package name: " + package_name)
deb_name = package_name + ".deb"

if os.path.exists(package_name):
    print("Aborting: build path \"" + package_name + "\" already exists (remove or create different version)")
    sys.exit(1)

if os.path.exists(deb_name):
    print("Aborting: The deb file to build (\""+ deb_name + "\") already exists")
    sys.exit(1)

# create basic package structure
os.makedirs(package_name + "/usr/lib")
os.makedirs(package_name + "/usr/include")
# copy library files
for module in modules:
    shutil.copy("../../../lib/sgpp/libsgpp" + module +".so", package_name + "/usr/lib")
# copy C++-headers
for module in modules:
    copy_headers.copy_headers(module, package_name)
# copy meta-files
os.makedirs(package_name + "/DEBIAN")
control_template = jinja2_env.get_template('control_template')
rendered = control_template.render(major_version=major_version, minor_version=minor_version, package_revision=package_revision, maintainer_name=maintainer_name, maintainer_email=maintainer_email)
with open(package_name + "/DEBIAN/control", "w") as f:
    # dpkg-deb requires newline at the end of the file
    f.write(rendered + "\n")

# finally, create the package
try:
    subprocess.check_call(["fakeroot", "dpkg-deb", "--build", package_name])
except Exception as e:
    print("Error: building the deb-file failed")
    print("reason: ")
    print(str(e))

# now create the python bindings package
package_name = "libsgpp-python_" + major_version + "." + minor_version + "-" + package_revision
print("python bindings package name: " + package_name)
deb_name = package_name + ".deb"

if os.path.exists(package_name):
    print("Aborting: build path \"" + package_name + "\" already exists (remove or create different version)")
    sys.exit(1)

if os.path.exists(deb_name):
    print("Aborting: The deb file to build (\""+ deb_name + "\") already exists")
    sys.exit(1)

os.makedirs(os.path.join(package_name, "usr/lib/python2.7/dist-packages/pysgpp/"))

# copy main shared library
shutil.copy("../../../lib/pysgpp/_pysgpp_swig.so", os.path.join(package_name, "usr/lib/python2.7/dist-packages/pysgpp/_pysgpp_swig.so"))
# copy library files
import copy_python
copy_python.copy_python(package_name)
# copy meta-files
os.makedirs(os.path.join(package_name, "DEBIAN"))
control_template = jinja2_env.get_template('control_python_template')
rendered = control_template.render(major_version=major_version, minor_version=minor_version, package_revision=package_revision, maintainer_name=maintainer_name, maintainer_email=maintainer_email)
with open(os.path.join(package_name, "DEBIAN/control"), "w") as f:
    # dpkg-deb requires newline at the end of the file
    f.write(rendered + "\n")

# finally, create the package
try:
    subprocess.check_call(["fakeroot", "dpkg-deb", "--build", package_name])
except Exception as e:
    print("Error: building the deb-file failed")
    print("reason: ")
    print(str(e))
