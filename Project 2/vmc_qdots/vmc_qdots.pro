TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp \
    split.cpp

HEADERS += \
    system.h \
    particle.h
LIBS += -LC:\Users\ramir\OneDrive\Documents\armadillo-7.900.1\examples\lib_win64 -llapack_win64_MT
LIBS += -LC:\Users\ramir\OneDrive\Documents\armadillo-7.900.1\examples\lib_win64 -lblas_win64_MT
INCLUDEPATH += C:\Users\ramir\OneDrive\Documents\armadillo-7.900.1\inc
DEPENDPATH += C:\Users\ramir\OneDrive\Documents\armadillo-7.900.1\inc
INCLUDEPATH += C:\Users\ramir\OneDrive\Documents\armadillo-7.900.1\armadillo_bits
DEPENDPATH += C:\Users\ramir\OneDrive\Documents\armadillo-7.900.1\armadillo_bits
