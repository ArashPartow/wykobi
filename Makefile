#
#  ***********************************************************************
#  *                                                                     *
#  * Wykobi Computational Geometry Library                               *
#  * Release Version 0.0.5                                               *
#  * http://www.wykobi.com                                               *
#  * Copyright (c) 2005-2017 Arash Partow, All Rights Reserved.          *
#  *                                                                     *
#  * The Wykobi computational library and its components are supplied    *
#  * under the terms of the General Wykobi License agreement. The        *
#  * contents of the Wykobi computational library and its components may *
#  * not be copied or disclosed except in accordance with the terms of   *
#  * that agreement.                                                     *
#  *                                                                     *
#  * URL: http://www.wykobi.com/license.html                             *
#  *                                                                     *
#  ***********************************************************************
#


COMPILER     = -c++
#COMPILER    = -clang
OPTIONS      = -pedantic-errors -Wall -Wextra -Werror -O3 -o
OPTIONS_LIBS = -pedantic-errors -Wall -Wextra -Werror -O3 -c
LINKER_OPT   = -L/usr/lib -lstdc++ -lm

CPP_SRC =

OBJECTS = $(CPP_SRC:.cpp=.o)

%.o: %.hpp %.cpp
	$(COMPILER) $(OPTIONS_LIBS) $*.cpp -o $@

all: $(OBJECTS) wykobi_build 

wykobi_build : wykobi_build.cpp wykobi.hpp wykobi_math.hpp $(OBJECTS)
	$(COMPILER) $(OPTIONS) wykobi_build wykobi_build.cpp $(OBJECTS) $(LINKER_OPT)

clean:
	rm -f core *.o *.bak *stackdump


#
# The End !
#
