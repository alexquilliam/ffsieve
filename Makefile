#partially based off of https://stackoverflow.com/questions/5178125/how-to-place-object-files-in-separate-subdirectory

CC := g++
TARGET := ffsieve

SRCDIR := src
INCDIR := $(shell find include -type d)
BUILDDIR := obj
SRCEXT := cpp

CFLAGS := -Wall -g -ftrapv -O2 -ftree-loop-vectorize -march=native
LIB := -lgmp -lgmpxx -l:libprofiler.so.0.4.8 -lpthread
INC := $(addprefix -I,$(INCDIR))

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

all: $(TARGET)

clean:
	@$(RM) -rf $(BUILDDIR)
	@rm $(TARGET)

-include $(OBJECTS:.o=.d)

#Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.d
	@cp -f $(BUILDDIR)/$*.d $(BUILDDIR)/$*.d.tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.o:|' < $(BUILDDIR)/$*.d.tmp > $(BUILDDIR)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.d.tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.d
	@rm -f $(BUILDDIR)/$*.d.tmp

.PHONY: all clean
