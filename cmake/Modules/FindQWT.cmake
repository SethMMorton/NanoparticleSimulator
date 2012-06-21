# - Try to find the Qwt includes and library
# which defines
#
# QWT_FOUND - system has Qwt
# QWT_INCLUDE_DIR - where to find qwt.h
# QWT_LIBRARIES - the libraries to link against to use Qwt
# QWT_LIBRARY - where to find the Qwt library (not for general use)
# QWT_6 - if setSamples is used instead of setData (QWT version 6)

# copyright (c) 2006 Thomas Moenicke thomas.moenicke@kdemail.net
#
# Redistribution and use is allowed according to the terms of the BSD license.

INCLUDE(CheckCXXSourceCompiles)

IF(NOT QT4_FOUND)
    INCLUDE(FindQt4)
ENDIF(NOT QT4_FOUND)

SET(QWT_FOUND FALSE)

IF(QT4_FOUND)

    SET(QWT_NAMES ${QWT_NAMES} qwt qwt-qt4)
    FIND_LIBRARY(QWT_LIBRARY
        NAMES ${QWT_NAMES}
        PATHS /usr/local/qwt/lib /usr/local/lib /usr/lib /usr/lib64
    )

    IF(QWT_LIBRARY)

        SET(QWT_LIBRARIES ${QWT_LIBRARY})
        SET(QWT_FOUND TRUE)

        IF(CYGWIN)
            IF(BUILD_SHARED_LIBS)
            # No need to define QWT_USE_DLL here, because it's default for
            # Cygwin.
            ELSE(BUILD_SHARED_LIBS)
            SET(QWT_DEFINITIONS -DQWT_STATIC)
            ENDIF(BUILD_SHARED_LIBS)
        ENDIF(CYGWIN)

    ENDIF(QWT_LIBRARY)

    FIND_PATH(QWT_INCLUDE_DIR qwt.h
    /usr/local/qwt/include
    /usr/local/include
    /usr/include/qwt-qt4
    /usr/include/qwt
    /usr/include
    ${QWT_LIBRARY}/Headers
    )

    # Check how to call set[Data|Samples]
    set(QWT_CXX_TEST_SOURCE 
"
#include<QtCore>
#include<qwt_plot_curve.h>
int main(int argc, char *argv[]) {

    QwtPlotCurve *temp = new QwtPlotCurve();
    double x[] = { 0.0, 1.0 };
    double y[] = { 0.0, 1.0 };
    temp->setSamples(x, y, 2);
    return 0;
}
")

    SET(CMAKE_REQUIRED_INCLUDES ${QWT_INCLUDE_DIR} ${QT_INCLUDES})
    SET(CMAKE_REQUIRED_LIBRARIES ${QWT_LIBRARIES} ${QT_LIBRARIES})
    UNSET(QWT_6_DETECTED CACHE)
    CHECK_CXX_SOURCE_COMPILES("${QWT_CXX_TEST_SOURCE}" QWT_6_DETECTED)
    IF(QWT_6_DETECTED)
        SET(QWT_6 TRUE)
    ENDIF(QWT_6_DETECTED)

ENDIF(QT4_FOUND)

IF(QWT_FOUND)
  IF(NOT QWT_FIND_QUIETLY)
    MESSAGE(STATUS "Found Qwt: ${QWT_LIBRARY}")
  ENDIF(NOT QWT_FIND_QUIETLY)
ELSE(QWT_FOUND)
  IF(QWT_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find Qwt library")
  ENDIF(QWT_FIND_REQUIRED)
ENDIF(QWT_FOUND)

MARK_AS_ADVANCED(QWT_INCLUDE_DIR QWT_LIBRARY, QWT_6)
