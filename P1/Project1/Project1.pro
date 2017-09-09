TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    Project1.cpp
LIBS += -larmadillo -llapack -lblas
