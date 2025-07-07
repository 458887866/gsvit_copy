#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.50])
AC_INIT(gsvit, 1.9.4, pklapetek@cmi.cz)
AC_CONFIG_SRCDIR([src3d/svit.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AM_PROG_CC_C_O

# Checks for libraries.
GLIB_REQUIRED_VERSION=2.0.0
AM_PATH_GLIB_2_0($GLIB_REQUIRED_VERSION,,[AC_MSG_ERROR(Test for GLib failed.)], glib)
GTK_REQUIRED_VERSION=2.0.0
AM_PATH_GTK_2_0($GTK_REQUIRED_VERSION,,[AC_MSG_ERROR(Test for GTk failed.)], gtk)


# Checks for Gwyddion
 AC_ARG_WITH([gwyddion],
    [  --with-gwyddion=PATH        prefix where gwyddion is installed [default=auto]])
 if test "x$with_gwyddion" != x; then
  GWYDDION_CFLAGS="-I$with_gwyddion/include/gwyddion/ -I$with_gwyddion/lib/gwyddion/include/ -Wall -I$with_gwyddion/lib64/gwyddion/include/"
  GWYDDION_LIBS="-L$with_gwyddion/lib -L$with_gwyddion/lib64 -lgwyprocess2 -lgwyddion2 -lgwymodule2 -lgwyapp2 -lgwydraw2 -lgwydgets2 -lgobject-2.0 -lm"
  echo "Gwyddion: path supplied by user used."
  GWYDDION_VERSION=$($with_gwyddion/bin/gwyddion --version) 
  if test "x$GWYDDION_VERSION" != x; then
     echo "Gwyddion succesfully found."
  else
     AC_MSG_ERROR(Test for Gwyddion failed.)   
  fi
 else
  PKG_CHECK_MODULES(GWYDDION, gwyddion)
  echo "Gwyddion: default path used."
 fi
GWYDDION_CFLAGS="$GWYDDION_CFLAGS -O3 -fopenmp -g -lm -I/usr/include/gtkglext-1.0 -I/usr/lib64/gtkglext-1.0/include/"
GWYDDION_LIBS="$GWYDDION_LIBS -lgthread-2.0 -lm -L/usr/lib64 -export-dynamic -lgtkglext-x11-1.0 -lgdkglext-x11-1.0 -lGLU -lGL -lXmu -lXt -lSM -lICE -lgtk-x11-2.0 -lpangox-1.0 -lX11 -lgdk-x11-2.0 -latk-1.0 -lgio-2.0 -lpangoft2-1.0 -lgdk_pixbuf-2.0 -lpangocairo-1.0 -lcairo -lpango-1.0 -lfreetype -lz -lfontconfig -lgobject-2.0 -lgmodule-2.0 -lglib-2.0 -lfftw3"
AC_SUBST(GWYDDION_CFLAGS)
AC_SUBST(GWYDDION_LIBS)

# Checks for CUDA
AC_ARG_WITH([cuda],
    [  --with-cuda=PATH        prefix where cuda is installed [default=auto]])
if test "x$with_cuda" != x; then
  CUDA_CFLAGS="-I$with_cuda/include -Wall"
  CUDA_LIBS="-L$with_cuda/lib -L$with_cuda/lib64 -lcudart"
  NVCC="$with_cuda/bin/nvcc"
  echo "CUDA: path supplied by user used."
else
  CUDA_CFLAGS="-I/usr/local/cuda/include -Wall"
  CUDA_LIBS="-L/usr/local/cuda/lib -L/usr/local/cuda/lib64 -lcudart"
  NVCC="nvcc"
  echo "CUDA: default path used."
fi

NVCCFLAGS="-O2 -use_fast_math -DUCUDA" 
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(NVCC)
AC_SUBST(NVCCFLAGS)

NVCCVERSION=$($NVCC --version)
if test "x$NVCCVERSION" = x; then
       echo "Cannot run NVIDIA compiler, configuring with no CUDA support."
       AC_DEFINE(USE_CUDA, 0, [Cuda compiler not found])
       AM_CONDITIONAL([USE_CUDA], [false])
       UCUDA_CFLAGS="-UUCUDA" 
else
       echo "NVCC version is $NVCCVERSION, configuring with CUDA support."
       AC_DEFINE(USE_CUDA, 1, [Cuda compiler installed.])
       AM_CONDITIONAL([USE_CUDA], [true])
       UCUDA_CFLAGS="-DUCUDA"

fi
AC_SUBST(USE_CUDA)
AC_SUBST(UCUDA_CFLAGS)


AC_CONFIG_FILES([Makefile
		 gsvit.spec:gsvit.spec.in
                 src2d/Makefile
                 src3d/Makefile
                 tests/Makefile
                 tests/test01/Makefile
                 tests/test02/Makefile
                 tests/test03/Makefile
                 tests/test04/Makefile
                 tests/test05/Makefile
                 tests/test06/Makefile
                 tests/selftests/Makefile
                 data/Makefile
                 data/spectra/Makefile
                 ])
AC_OUTPUT
