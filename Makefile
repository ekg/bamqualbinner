HEADERS = 
SOURCES = bamqualbinner.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BINS = bamqualbinner

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64
INCLUDES = -I$(BAMTOOLS_ROOT)/src
LDFLAGS =
LIBS = -lz -lm -L./ -lbamtools

BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib

# profiling

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

# libraries

# builds bamtools static lib, and copies into root
libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)
	cp bamtools/lib/libbamtools.a ./

# bamqualbinner build

%.o: %.cpp %.h
	$(CXX) -c -o $@ $(*F).cpp $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(BINS): $(BIN_SOURCES) $(OBJECTS) $(SOURCES) $(HEADERS) libbamtools.a
	$(CXX) $(OBJECTS) -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(BINS) $(OBJECTS)
	cd bamtools/build && $(MAKE) clean
	rm libbamtools.a

clean-bamqualbinner:
	rm -f $(BINS) $(OBJECTS)

.PHONY: clean all
