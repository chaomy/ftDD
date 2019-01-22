

# compiler 
CXX       =  g++
CPPFLAGS +=  -Wall -g -O3 -std=c++11

# include and link  
CINCLUDE  =  -I../include 
CLINK 	  =   
CDLINK    =  -L${CLINK} -lm -larmadillo -framework Accelerate

# source 
SRCS = ftDDMain.cpp \
	   ftDDload.cpp \
	   ftDDinit.cpp \
	   ftDDcmaes.cpp \
	   ftDDerr.cpp  \
	   ftDDerrFix.cpp \
	   ftDDopt.cpp  \
	   ftDDParam.cpp \
	   ftDDpreprocess.cpp \
	   ftDDout.cpp  	    

# target
MAKETARGET 	= ftDD.exe

OBJECTS := $(subst .cpp,.o,${SRCS})

CC = ${CXX}

%.o: %.cpp
	${CC} -c  $<  ${CINCLUDE} ${CPPFLAGS}

${MAKETARGET}:$(OBJECTS)

all: ${MAKETARGET}
	${CC}  -o  ${MAKETARGET}  ${OBJECTS}  ${CDLINK}  ${CPPFLAGS}

clean:
	rm -f *.o