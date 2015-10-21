// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _GMSH_CONFIG_H_
#define _GMSH_CONFIG_H_

/* #undef HAVE_3M */
/* #undef HAVE_64BIT_SIZE_T */
/* #undef HAVE_ACIS */
#define HAVE_ANN
#define HAVE_BAMG
#define HAVE_BFGS
#define HAVE_BLAS
#define HAVE_BLOSSOM
#define HAVE_CHACO
#define HAVE_DLOPEN
#define HAVE_DINTEGRATION
#define HAVE_FLTK
#define HAVE_FL_TREE
/* #undef HAVE_FOURIER_MODEL */
#define HAVE_GMM
#define HAVE_GMP
#define HAVE_KBIPACK
#define HAVE_LAPACK
/* #undef HAVE_LIBCGNS */
#define HAVE_LIBJPEG
#define HAVE_LIBPNG
#define HAVE_LIBZ
#define HAVE_MATHEX
#define HAVE_MED
#define HAVE_MESH
#define HAVE_METIS
#define HAVE_MMG3D
#define HAVE_MPEG_ENCODE
/* #undef HAVE_MPI */
#define HAVE_NATIVE_FILE_CHOOSER
#define HAVE_NETGEN
/* #undef HAVE_NO_INTPTR_T */
/* #undef HAVE_NO_SOCKLEN_T */
/* #undef HAVE_NO_STDINT_H */
/* #undef HAVE_NO_VSNPRINTF */
#define HAVE_OCC
#define HAVE_ONELAB
#define HAVE_OPENGL
/* #undef HAVE_OSMESA */
#define HAVE_PARSER
#define HAVE_PETSC
#define HAVE_PLUGINS
#define HAVE_POST
/* #undef HAVE_QT */
#define HAVE_RTREE
#define HAVE_SALOME
#define HAVE_SLEPC
#define HAVE_SOLVER
/* #undef HAVE_TAUCS */
#define HAVE_TETGEN
#define HAVE_VORO3D

#define GMSH_CONFIG_OPTIONS " Ann Bamg Bfgs Blas(VecLib) Blossom Chaco DIntegration Dlopen FlTree Fltk GMP Gmm Jpeg(Fltk) Kbipack Lapack(VecLib) MathEx Med Mesh Metis Mmg3d Mpeg NativeFileChooser Netgen OneLab OpenCascade OpenGL OptHom PETSc Parser Plugins Png(Fltk) Post RTree SLEPc Salome Solver Tetgen(1.5) Voro3D Zlib"



#endif
