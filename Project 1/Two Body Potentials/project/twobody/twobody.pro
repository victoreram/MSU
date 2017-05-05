TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    ../../Two_body_potentials/elements.cpp \
    ../../Two_body_potentials/Coulomb_Functions.cpp

HEADERS += \
    ../../Two_body_potentials/Coulomb_Functions.hpp
LIBS += -LC:/Users/ramir/OneDrive/Documents/armadillo-7.700.0/armadillo-7.700.0/examples/lib_win64 -lblas_win64_MT
LIBS += -LC:/Users/ramir/OneDrive/Documents/armadillo-7.700.0/armadillo-7.700.0/examples/lib_win64 -llapack_win64_MT

INCLUDEPATH += C:/Users/ramir/OneDrive/Documents/armadillo-7.700.0/armadillo-7.700.0/include
DEPENDPATH += C:/Users/ramir/OneDrive/Documents/armadillo-7.700.0/armadillo-7.700.0/include
