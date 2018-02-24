TARGET = ffsieve
LIBS = -lgmp
CC = g++
CFLAGS = -Wall -g

.PHONY: default all clean

default: $(TARGET)
all: default

INC_DIRS := $(shell find include -type d)
	C_INC_FLAGS += $(addprefix -I,$(INC_DIRS))

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(C_INC_FLAGS) -c $< -o $@

OBJECTS = $(patsubst %.cpp, %.o, $(shell find . -name '*.cpp'))
HEADERS = $(shell find . -name '*.hpp')

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	find . -type f -name '*.o' -exec rm {} +
	find . -type f -name $(TARGET) -exec rm {} +
