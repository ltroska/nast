from distutils.core import setup
from distutils.extension import Extension

import os

os.environ["CFLAGS"] += ' -I/home/ltroska/.buildozer/android/platform/android-ndk-r17c/sources/cxx-stl/stlport/stlport'
os.environ["CXXFLAGS"] += ' -I/home/ltroska/.buildozer/android/platform/android-ndk-r17c/sources/cxx-stl/stlport/stlport'
os.environ["LDFLAGS"] += ' -L/home/ltroska/.buildozer/android/platform/android-ndk-r17c/sources/cxx-stl/stlport/libs/armeabi-v7a -lstlport_shared'

ext_modules = [Extension("nast", ["nast.cpp"], language="c++")]

setup(name='nast', ext_modules = ext_modules)
