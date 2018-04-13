// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

package sgpp;

public class LoadJSGPPLib {
  public static void loadJSGPPLib() {
    // load jsgpp library
    try {
      System.loadLibrary("jsgpp");
    } catch (UnsatisfiedLinkError e) {
      // under win32, loadLibrary looks after <NAME>.dll,
      // so we have to add our "lib" prefix
      System.loadLibrary("libjsgpp");
    }
  }
}
