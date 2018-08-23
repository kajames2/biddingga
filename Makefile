CC := clang++-6.0
SRCDIR := src
BUILDDIR := build
TESTDIR := test
TARGET := bidding_ga

SRCEXT := cc
SOURCES := $(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
TESTS := $(shell find $(TESTDIR) -type f -name "*.$(SRCEXT)")
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TESTOBJECTS := $(patsubst $(TESTDIR)/%,$(BUILDDIR)/%,$(TESTS:.$(SRCEXT)=.o))
CFLAGS := -g -std=c++17 -Wl,-rpath,./lib/ -O3 -flto #-fopenmp
LIB := -lpthread -lbidding -lgenericga -Llib/
INC := -I include -I lib

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo "$(CC) $(CFLAGS) $(INC) $^ $(LIB) $(TARGET).$(SRCEXT) -o bin/$(TARGET)"; $(CC) $(CFLAGS) $(INC) $^ $(LIB) $(TARGET).$(SRCEXT) -o bin/$(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(@D)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Tests
tests_main: $(TESTOBJECTS) $(OBJECTS)
	$(CC) $(CFLAGS) tests_main.cc $(INC)  $(TESTOBJECTS) $(OBJECTS) -o bin/tests_main $(LIB)

$(BUILDDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	@mkdir -p $(@D)
	@echo " Making $(@D)"
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

.PHONY: clean
