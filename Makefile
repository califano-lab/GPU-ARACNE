CC := /usr/local/cuda/bin/nvcc # This is the main compiler

SRCDIR := src
BUILDDIR := build
TARGET := bin/runner  
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -Xcompiler "-std=c++0x" -g  
LIB := 
INC := -I include -arch=sm_30

$(TARGET): $(OBJECTS)
	  @echo " Linking..."
	    @echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	  @mkdir -p $(BUILDDIR)
	    @echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

clean:
	  @echo " Cleaning..."; 
	    @echo " $(RM) -r $(BUILDDIR) $(TARGET) *~"; $(RM) -r $(BUILDDIR) $(TARGET) *~

.PHONY: clean