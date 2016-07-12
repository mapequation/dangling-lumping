# Various flags
CXX  = clang++
LINK = $(CXX)
# CXXFLAGS = -I -Wall -g
CXXFLAGS = -std=c++11 -Wall -O3
LFLAGS = -lm

TARGET  = dangling-lumping

HEADER  = dangling-lumping.h
FILES = dangling-lumping.cc

OBJECTS = $(FILES:.cc=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $^ $(LFLAGS) -o $@

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile




