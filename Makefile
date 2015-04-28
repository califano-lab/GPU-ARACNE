CC := /usr/local/cuda/bin/nvcc # This is the main compiler

SRCDIR := src
BUILDDIR := build
TARGET := bin/runner  
SRCEXT := cu
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g #-Xcompiler "-std=c++0x" -g  
LIB := 
INC := -I src -I include -arch=sm_35

$(TARGET): $(OBJECTS)
	  @echo " Linking..."
	    @echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	  @mkdir -p $(BUILDDIR)
	    @echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cu  $(INC) $(LIB) -o bin/tester
jTester:
	$(CC) $(CFLAGS) test/jingTester.cu  $(INC) $(LIB) -o bin/jingTester
unitTester: 
	$(CC) $(CFLAGS) test/unitTest_DPI.cu  $(INC) $(LIB) -o bin/unitTester
clean:
	  @echo " Cleaning..."; 
	    @echo " $(RM) -r $(BUILDDIR) $(TARGET) src/*~"; $(RM) -r $(BUILDDIR) $(TARGET) src/*~ *~

.PHONY: clean
