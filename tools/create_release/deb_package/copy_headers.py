#!/usr/bin/python3
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import re
import os
import shutil

def copy_headers(module, package_name):
    module_src_folder = "../../../" + module + "/src/"
    # print(module_src_folder)
    # print("package_name: " + package_name)
    package_realpath = os.path.realpath(package_name)
    # print(os.path.abspath(package_name))
    for root, dirs, files in os.walk(module_src_folder):
        for f in files:
            if not f.endswith(".hpp"):
                continue
            origin_path = os.path.join(root, f)
            # print("origin:" + origin_path)
            target_suffix = re.search(r"\.\./\.\./\.\./" + module + r"/src/" + r"(.*)", origin_path).group(1)
            target_path = os.path.join(package_realpath, "usr/include", target_suffix)
            target_dir = os.path.dirname(target_path)
            print("copy: " + origin_path + " -> " + target_path)
            os.makedirs(target_dir, exist_ok=True)
            shutil.copy(origin_path, target_path)
