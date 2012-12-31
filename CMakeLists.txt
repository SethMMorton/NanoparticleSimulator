# CMAKE file for the NPSim project

####################
# SET UP ENVIRONMENT
####################

PROJECT(NPSim CXX C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

##################################
# DEFINE WHAT IS NEEDED TO COMPILE
##################################

# We need QT
FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtSvg)
INCLUDE(${QT_USE_FILE})
# We also need QWT
FIND_PACKAGE(QWT REQUIRED)

##############################################
# DEFINE THE FILES AND FOLDERS IN THIS PACKAGE
##############################################

# The source and include directories
SET(SRC src)
SET(INC include)
SET(LIB lib)

# Directories containing include files
INCLUDE_DIRECTORIES(${INC}
                    ${QWT_INCLUDE_DIR})

# The linker should look for libraries in lib
LINK_DIRECTORIES(${LIB})

############################################
# DEFINE HOW TO COMPILE THE NPSIM EXECUTABLE
############################################

# The name of the NPSim executable
SET(NPSIMEXE NPSim${BINARY_POSTFIX})

# Add the headers and source files for the main program
SET(NPSim_HEADERS ${INC}/mainwindow.hpp
                  ${INC}/nanoparticle.hpp
                  #${INC}/QScienceSpinBox.hpp
                  #${INC}/plot.hpp
                  #${INC}/inputfunction.hpp
                  #${INC}/function.hpp
                  #${INC}/xypoints.hpp
                  #${INC}/valuebox.hpp
)

SET(NPSim_SRC ${SRC}/main.cpp
              ${SRC}/nanoparticle.cpp
              ${SRC}/mainwindow.cpp
              #${SRC}/QScienceSpinBox.cpp
              #${SRC}/plot.cpp
              #${SRC}/inputfunction.cpp
              #${SRC}/function.cpp
              #${SRC}/xypoints.cpp
              #${SRC}/valuebox.cpp
              #${SRC}/defaultfunctions.cpp
)

# Define the files that need to be wrapped with moc
SET(NPSim_MOC_HEADERS #${INC}/mainwindow.hpp
                      #${INC}/QScienceSpinBox.hpp
                      #${INC}/plot.hpp
                      #${INC}/inputfunction.hpp
                      #${INC}/function.hpp
                      #${INC}/valuebox.hpp
)


# Wrap with moc
QT4_WRAP_CPP(NPSim_MOC ${NPSim_MOC_HEADERS})

# Add a definition to differentiate versions of QWT
IF (QWT_6_DETECTED)
    ADD_DEFINITIONS(-DQWT_6_DETECTED)
ENDIF (QWT_6_DETECTED)

# Define the files to be used to make the executable
ADD_EXECUTABLE(${NPSIMEXE}
               ${NPSim_SRC}
               ${NPSim_HEADERS}
               ${NPSim_MOC})

# Link executable to Qt
TARGET_LINK_LIBRARIES(${NPSIMEXE}
                      ${NPSOLVELIB}
                      ${QT_LIBRARIES}
                      ${QWT_LIBRARIES})
