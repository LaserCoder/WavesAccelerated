cmake_minimum_required(VERSION 3.5)

project(NUFOSAleph VERSION 0.1 LANGUAGES CUDA CXX)
enable_language(CUDA)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Ofast -DEIGEN_USE_BLAS -lblas -llapack -llapacke -lfftw3 -lm")

set(CMAKE_PREFIX_PATH "/home/m3litaha/Qt/6.5.3/gcc_64/lib/cmake")
find_package(QT NAMES Qt6 Qt5  COMPONENTS Widgets PrintSupport REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets PrintSupport)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
        qcustomplot.cpp
        qcustomplot.h
        kernel.cu
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(NUFOSAleph
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET NUFOSAleph APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(NUFOSAleph SHARED
            ${PROJECT_SOURCES}
            qcustomplot.h
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(NUFOSAleph
            ${PROJECT_SOURCES}
            qcustomplot.h
        )
    endif()
endif()
set_target_properties(NUFOSAleph PROPERTIES CUDA_ARCHITECTURES 86)
set_target_properties(NUFOSAleph PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(NUFOSAleph PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::PrintSupport Eigen3::Eigen)
target_link_libraries(NUFOSAleph PRIVATE fftw3)
# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.
if(${QT_VERSION} VERSION_LESS 6.1.0)
  set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.NUFOSAleph)
endif()
set_target_properties(NUFOSAleph PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS NUFOSAleph
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)



if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(NUFOSAleph)
endif()
