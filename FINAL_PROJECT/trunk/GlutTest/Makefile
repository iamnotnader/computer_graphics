# 
# List of source files
#
MESHVIEW_SRCS= GlutTest.cpp R3Scene.cpp imageloader.cpp vec3f.cpp smoke.cpp  
MESHVIEW_OBJS=$(MESHVIEW_SRCS:.cpp=.o)


#
# Compile and link options
#

CXX=g++
CXXFLAGS=-Wall -I. -g -DUSE_JPEG


#
# OpenGL libraries
#
UNAME := $(shell uname)
ifneq (,$(findstring Darwin,$(UNAME)))
	GLLIBS = -framework GLUT -framework OpenGL
else
  ifneq (,$(findstring CYGWIN,$(UNAME)))
	GLLIBS = -lopengl32 -lglu32 -lglut32
  else
	GLLIBS = -lglut -lGLU -lGL
  endif
endif



#
# Compile command
#

%.o: %.cpp R3Mesh.h
	    $(CXX) $(CXXFLAGS) -c $< -o $@

#
# GNU Make: targets that don't build files
#

.PHONY: all clean distclean

#
# Rules encoding targets and dependencies.  By default, the first of
# these is built, but you can also build any individual target by
# passing it to make - e.g., "make imgpro" or "make clean"
#
# Notice that many of the dependencies are implicit (e.g. a .o depends
# on its corresponding .cpp), as are many of the compilation rules.
#

LIBS=jpeg/libjpeg.a R3/libR3.a R2/libR2.a

all: starfox

R3/libR3.a: 
	    $(MAKE) -C R3

R2/libR2.a: 
	    $(MAKE) -C R2

jpeg/libjpeg.a: 
	    $(MAKE) -C jpeg

starfox: $(LIBS) $(MESHVIEW_OBJS)
	rm -f $@
	$(CXX) $(CXXFLAGS) $^ -lm -o $@ $(LIBS) $(GLLIBS)

clean: 
	rm -f *.o meshpro starfox
	$(MAKE) -C R3 clean
	$(MAKE) -C R2 clean
	$(MAKE) -C jpeg clean
