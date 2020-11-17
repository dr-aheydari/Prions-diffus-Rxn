# to run with OMP on macOS
macx:
{
QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include
}
macx:
{
QMAKE_LFLAGS += -lomp
}

macx:
{
LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib
}


# mains
SOURCES += \
#example_mains/QuadTreeSolverNodeBasedPoissonJump_example_main.cpp\
#example_mains/OcTreeSolverFacesPoissonJump_example_main.cpp\
#example_mains/QuadTreeSolverCellCenteredPoissonJump_example_main.cpp\
#example_mains/QuadTreeSolverFacesPoissonJump_example_main.cpp\
#example_mains/OcTreeSolverCellCenteredPoissonJump_example_main.cpp\
#example_mains/Advection_SemiLagrangian_example_main.cpp\
#example_mains/Advection_Augmented_LS_example_main.cpp\
#example_mains/main_test_cubic_interpolation.cpp\
#mains_ElectroHydro/main_ElectroHydro.cpp\
#mains_ElectroHydro/main_ElectroHydro.cpp\
#mains_ElectroHydro/main_ElectroHydro_SolidSphere.cpp\
#mains_ElectroHydro/main_ElectroHydro_backup.cpp\
#mains_ElectroHydro/main_ElectroHydro_3DSolidSphere.cpp\
#mains_ElectroHydro/main_ElectroHydro_3D_2drops.cpp\
#mains_ElectroHydro/main_ElectroHydro_3D_SolidSphere.cpp\
#mains_ElectroHydro/main_ElectroHydro_test_Godunov.cpp\
#mains_NavierStokes/Diphasic/main_test_extension_velocity.cpp\
# mains_NavierStokes/Diphasic/main_NS_diphasic_oscillating_bubble.cpp\
#mains_NavierStokes/Diphasic/test_OpenMP.cpp\
#mains_NavierStokes/Diphasic/main_NS_diphasic_oscillating_bubble_3D.cpp\
#mains_NavierStokes/mains_for_maxime/main_NS_Karman_vortex_street.cpp\
#mains_NavierStokes/Diphasic/main_NS_diphasic_parrasite_currents.cpp\
#mains_NavierStokes/Diphasic/main_NS_diphasic_parrasite_current_3D.cpp\
#mains_NavierStokes/Diphasic/main_NS_diphasic_analityc.cpp \
mains_prions/prion_cell_splitting.cpp
#mains_NavierStokes/Diphasic/main_NS_diphasic_rising_bubbles_3D.cpp\
#mains_NavierStokes/Diphasic/main_NS_diphasic_bubbler_3D.cpp\
#example_mains/main_test_cubic_interpolation_3D.cpp\
#mains_NavierStokes/Diphasic/main_NS_diphasic_rising_bubble.cpp\


TEMPLATE = app
TARGET = Prions_Project

# MAC MODIFIED FOR OPENMP
### OMP,PETSC STUFF: only for release...
message($${CONFIG})
!contains(CONFIG, "gcc debug"){
message( "This is a message" )
#QMAKE_CXXFLAGS += -fopenmp
#QMAKE_LFLAGS += -fopenmp # -lvoro++
DEFINES += CASL_THROWS\
CASL_PETSC\
#CASL_OPENMP\
#LIBS += -fopenmp
}
DEFINES += CASL_THROWS\

QMAKE_CXXFLAGS +=  -w -O3
QMAKE_LFLAGS +=  -w -O3
#LIBS += -fopenmp

#LIBS += -lvoro++

###
PETSC_DIR =  /Users/aliheydari/Documents/Softwares/petsc-3.10.0
PETSC_ARCH = arch-darwin-cxx-debug
PETSC_INCLUDES = $${PETSC_DIR}/include $${PETSC_DIR}/$${PETSC_ARCH}/include
PETSC_LIBS = -L$${PETSC_DIR}/$${PETSC_ARCH}/lib -Wl,-rpath,$${PETSC_DIR}/$${PETSC_ARCH}/lib -lpetsc -lmpich # -L/usr/X11R6/lib  -lHYPRE -lpthread -lflapack -lfblas  -ldl -lmpichf90 -lpthread -lquadmath  -lm  -lstdc++  -lopa -lmpl -lpetsc
#LIBS += -L/$${PETSC_DIR}/$${PETSC_ARCH}/lib -lpetsc -L/usr/X11R6/lib -lHYPRE -lpthread -lflapack -lfblas  -ldl -lmpichf90 -lpthread -lquadmath -l -lm -lstdc++  -lmpich -lopa -lmpl  -ldl


#LIBS += $${PETSC_DIR}/$${PETSC_ARCH}/lib

## this one was obtain with make petscdirc petsc arch getlinklibs

LIBS+=   -Wl,-rpath,/$${PETSC_DIR}/$${PETSC_ARCH}/lib # -L/home/maxime/soft/petsc-3.5.2/arch-linux2-cxx-opt/lib
# -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8
#-Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu #-L/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl #-Wl,-rpath,/home/maxime/soft/petsc-3.5.2/arch-linux2-cxx-opt/lib #-lmpich -lopa -lmpl -lrt -lpthread -lgcc_s -ldl


LIBS += $${PETSC_LIBS}

CONFIGURE_ARGS += -vt

DEFINES += CASL_MAX_IS_USING_IT\
#CASL_THROWS \
CASL_TRUE_BC_HODGE\
#CASL_CELL_TRUNCATION\
CASL_TRUE_BC_VELO\
#CASL_UNSTABLE_DIV





#CASL_PATH = /home/maxime/Documents/Codes/CASL_Library
CASL_PATH = /home/maxime/Documents/Codes/CASL_Library
VORO_PATH = /Users/aliheydari/Documents/Softwares/voro++-0.4.6/src

#INCLUDEPATH += $${CASL_PATH}
DEPENDPATH += $${CASL_PATH}
INCLUDEPATH += $${PETSC_INCLUDES}
INCLUDEPATH += $${VORO_PATH}



# mains
SOURCES += lib/amr/OcTreeCellBasedLevelSet.cpp \
    lib/amr/OctNgbdNodesOfNode.cpp \
    lib/geometry/Cube3.cpp \
    lib/fastmarching/OcTreeFastMarching.cpp \
    lib/tools/OcTreeToolbox.cpp \
    lib/geometry/MarchingSquares.cpp \
    lib/algebra/PetscLinearSolver.cpp \
    lib/solvers/QuadTreeSolverNodeBasedPoissonJump.cpp \
#    lib/solvers/OcTreeSolverFacesPoissonJump_voro++.cpp \
    lib/amr/OcTreeCellBasedMAC.cpp \
    lib/amr/OcTreeCellBasedMACLevelSet.cpp \
#    lib/solvers/OcTreeSolverFacesPoisson_voro++.cpp \
#    lib/solvers/OcTreeSolverDiphasicFlow.cpp \
    lib/amr/OcTreeAdvection.cpp \
#    lib/solvers/OcTreeSolverNavierStokes.cpp
#    lib/amr/QuadTreeAugmentedLevelSet.cpp


# StopWatch
HEADERS += \
    lib/tools/StopWatch.h \
    lib/algebra/Cholesky.h \
    lib/amr/OcTreeLevelSet.h \
    lib/amr/OctNgbdNodesOfNode.h \
    lib/geometry/Cube3.h \
    lib/fastmarching/OcTreeFastMarching.h \
    lib/tools/LinearInterpolationOnOcTree.h \
    lib/tools/QuadraticInterpolationOnOcTree.h \
    lib/tools/OcTreeToolbox.h \
#    lib/solvers/OcTreeSolverNavierStokes.h \
#    lib/solvers/OcTreeSolverFacesPoisson_voro++.h \
    lib/amr/OcTreeCellBasedMAC.h \
    lib/amr/OcTreeCellBasedLevelSet.h \
    lib/amr/OcTreeAdvection.h \
    lib/solvers/OcTreeSolverCellCenteredPoisson.h \
    lib/amr/OcTreeCellBasedMACLevelSet.h \
    lib/geometry/MarchingSquares.h \
    lib/algebra/PetscLinearSolver.h \
    lib/solvers/QuadTreeSolverNodeBasedPoissonJump.h \
#   lib/solvers/OcTreeSolverFacesPoissonJump_voro++.h \
#    lib/solvers/OcTreeSolverDiphasicFlow.h \
#    lib/solvers/QuadTreeSolverActiveFluid.h
#    lib/amr/QuadTreeAugmentedLevelSet.h


SOURCES += \
   lib/tools/StopWatch.cpp


# Arrays
HEADERS += \
    lib/arrays/ArrayV.h \
    lib/arrays/ArrayV2D.h \
    lib/arrays/ArrayV3D.h
SOURCES += \
    lib/arrays/ArrayV.cpp \
    lib/arrays/ArrayV2D.cpp \
    lib/arrays/ArrayV3D.cpp


# Geometry
HEADERS += \
    lib/geometry/Point2.h \
    lib/geometry/Point3.h \
    lib/geometry/Simplex2.h \
    lib/geometry/Cube2.h
SOURCES += \
    lib/geometry/Point2.cpp \
    lib/geometry/Point3.cpp \
    lib/geometry/Simplex2.cpp \
    lib/geometry/Cube2.cpp

# AMR
HEADERS += \
    lib/amr/QuadNode.h \
    lib/amr/QuadCell.h \
    lib/amr/QuadTree.h \
    lib/amr/QuadTreeCellBased.h \
    lib/amr/QuadTreeCellBasedMAC.h \
    lib/amr/QuadNgbdNodesOfNode.h \
    lib/amr/QuadTreeLevelSet.h \
    lib/amr/OcTreeLevelSet.h \
    lib/amr/QuadTreeCellBasedLevelSet.h \
    lib/amr/QuadTreeCellBasedMACLevelSet.h \
    lib/amr/QuadTreeAdvection.h \
    lib/amr/OcTree.h \
    lib/amr/OctCell.h \
    lib/amr/OcTreeCellBased.h \

SOURCES += \
    lib/amr/QuadCell.cpp \
    lib/amr/QuadTree.cpp \
    lib/amr/QuadTreeCellBased.cpp \
    lib/amr/QuadTreeCellBasedMAC.cpp \
    lib/amr/QuadNgbdNodesOfNode.cpp \
    lib/amr/QuadTreeLevelSet.cpp \
    lib/amr/OcTreeLevelSet.cpp \
    lib/amr/QuadTreeCellBasedLevelSet.cpp \
    lib/amr/QuadTreeCellBasedMACLevelSet.cpp \
    lib/amr/QuadTreeAdvection.cpp \
    lib/amr/OcTree.cpp \
    lib/amr/OctCell.cpp \
    lib/amr/OcTreeCellBased.cpp \

# Tools
HEADERS += \
    lib/tools/CASL_math.h \
    lib/tools/CASL_types.h \
    lib/tools/QuadTreeToolbox.h \
    lib/tools/LinearInterpolationOnQuadTree.h \
  lib/tools/QuadraticInterpolationOnQuadTree.h
SOURCES += \
    lib/tools/CASL_math.cpp \
    lib/tools/QuadTreeToolbox.cpp \
    lib/tools/LinearInterpolationOnQuadTree.cpp \
    lib/tools/QuadraticInterpolationOnQuadTree.cpp

# Algebra
HEADERS += \
    lib/algebra/Matrix.h \
    lib/algebra/MatrixFull.h \
    lib/algebra/SparseMatrix.h \
    lib/algebra/SparseMatrix_CRS.h \
    lib/algebra/BCGSTAB.h \
    lib/algebra/CG.h \
    lib/algebra/Cholesky.h\
    lib/algebra/Multigrid.h \
    lib/algebra/LinearSolver.h
SOURCES += \
    lib/algebra/Matrix.cpp \
    lib/algebra/MatrixFull.cpp \
    lib/algebra/SparseMatrix.cpp \
    lib/algebra/SparseMatrix_CRS.cpp \
    lib/algebra/BCGSTAB.cpp \
    lib/algebra/CG.cpp \
    lib/algebra/Cholesky.cpp\
    lib/algebra/Multigrid.cpp \
    lib/algebra/LinearSolver.cpp

# Solvers
HEADERS += \
    lib/solvers/QuadTreeSolverCellCenteredPoisson.h \
    lib/solvers/QuadTreeSolverCellCenteredPoissonJump.h \
    lib/solvers/OcTreeSolverCellCenteredPoisson.h \
    lib/solvers/OcTreeSolverCellCenteredPoissonJump.h \
    lib/solvers/QuadTreeSolverNodeBasedPoisson.h \
    lib/solvers/OcTreeSolverNodeBasedPoissonJump.h \
    lib/solvers/OcTreeSolverNodeBasedPoisson.h \
    lib/solvers/QuadTreeSolverNavierStokes.h \
    lib/solvers/QuadTreeSolverDiphasicFlow.h \
     lib/solvers/QuadTreeSolverFacesPoisson.h \
     lib/solvers/QuadTreeSolverFacesPoissonJump.h\
lib/solvers/QuadTreeSolverElectroHydro.h

SOURCES += \
    lib/solvers/QuadTreeSolverCellCenteredPoisson.cpp \
    lib/solvers/QuadTreeSolverCellCenteredPoissonJump.cpp \
lib/solvers/OcTreeSolverCellCenteredPoisson.cpp \
lib/solvers/OcTreeSolverCellCenteredPoissonJump.cpp \
    lib/solvers/QuadTreeSolverNodeBasedPoisson.cpp \
lib/solvers/OcTreeSolverNodeBasedPoissonJump.cpp \
lib/solvers/OcTreeSolverNodeBasedPoisson.cpp \
     lib/solvers/QuadTreeSolverNavierStokes.cpp \
lib/solvers/QuadTreeSolverDiphasicFlow.cpp \
     lib/solvers/QuadTreeSolverFacesPoisson.cpp \
     lib/solvers/QuadTreeSolverFacesPoissonJump.cpp\
lib/solvers/QuadTreeSolverElectroHydro.cpp


# Interpolations
HEADERS += \
    lib/tools/VInterpolationOnQuadTreeCellBasedMAC.h \
    lib/tools/UInterpolationOnQuadTreeCellBasedMAC.h\
#    lib/tools/QuadraticInterpolationOnQuadTree.cpp
    lib/tools/LinearInterpolationOnOcTree.h \
    lib/tools/QuadraticInterpolationOnOcTree.h
SOURCES += \
    lib/tools/VInterpolationOnQuadTreeCellBasedMAC.cpp \
    lib/tools/UInterpolationOnQuadTreeCellBasedMAC.cpp \
#    lib/tools/QuadraticInterpolationOnQuadTree.h
    lib/tools/LinearInterpolationOnOcTree.cpp \
    lib/tools/QuadraticInterpolationOnOcTree.cpp

HEADERS += \
    lib/geometry/Voronoi2D.h

SOURCES += \
    lib/geometry/Voronoi2D.cpp

#FMM
HEADERS += \
   lib/fastmarching/QuadTreeFastMarching.h

SOURCES += \
   lib/fastmarching/QuadTreeFastMarching.cpp

