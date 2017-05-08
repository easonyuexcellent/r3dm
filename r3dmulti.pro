QT += core
QT -= gui

CONFIG += c++11

TARGET = r3dmulti
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    basisfctn.cpp \
    config_class.cpp \
    gradient_class.cpp \
    hot_spots.cpp \
    matrix2d.cpp \
    OPTIMIZATION.cpp \
    register3d.cpp \
    splcoord.cpp \
    volume.cpp \
    opt_hot_spots.cpp \
    called_grd_max.cpp \
    called_grd_diff.cpp \
    node.cpp \
    called_grd_build.cpp \
    called_crd_upsp.cpp

DISTFILES += \
    r3dmulti.pro.user

HEADERS += \
    basisfctn.h \
    config_class.h \
    gradient_class.h \
    hot_spots.h \
    matrix2d.h \
    OPTIMIZATION.h \
    register3d.h \
    splcoord.h \
    volume.h \
    head.h \
    opt_hot_spots.h \
    called_grd_max.h \
    called_grd_diff.h \
    node.h \
    called_grd_build.h \
    called_crd_upsp.h
