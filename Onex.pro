TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += src/main.cpp \
    src/Algo.cpp \
    src/deque.cpp \
    src/Grouping.cpp \
    src/Groups.cpp \
    src/online.cpp \
    src/OnlineSession.cpp \
    src/TimeSeries.cpp \
    src/trillionDTW.cpp \
    src/util.cpp
    

HEADERS += \
    src/Algo.h \
    src/deque.h \
    src/Grouping.h \
    src/Groups.h \
    src/online.h \
    src/OnlineSession.h \
    src/TimeSeries.h \
    src/trillionDTW.h \
    src/util.h
