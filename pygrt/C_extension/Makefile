
SRC_DIR := src
INC_DIR := include
BUILD_DIR := build
BIN_DIR := bin
LIB_DIR := lib
LNAME := grt
LIB_NAME := lib$(LNAME).so

VERSION := $(shell cat version )
CURRENT_DIR := $(shell pwd)

CC := gcc
FOMPFLAGS := -fopenmp 
#  -fdump-tree-all  -g -ffast-math -O3 -fno-associative-math  -march=native -mtune=native 
CFLAGS :=  -O3 -g \
          -fPIC -Wall -Wextra -I$(INC_DIR) \
		  -lm -lfftw3  -lfftw3f -DGRT_VERSION="\"$(VERSION)\""
LDFLAGS := -shared -lm
RPATH_FLAGS :=  -L$(LIB_DIR) -l$(LNAME)

INCS := $(wildcard $(INC_DIR)/*.h)
SRCS := $(wildcard $(SRC_DIR)/*.c)
SRCS_NOMAIN := $(SRCS)
# remove main functions
SRCS_NOMAIN := $(filter-out $(SRC_DIR)/grt_main.c, $(SRCS_NOMAIN))
SRCS_NOMAIN := $(filter-out $(SRC_DIR)/grt_syn.c, $(SRCS_NOMAIN))
SRCS_NOMAIN := $(filter-out $(SRC_DIR)/grt_travt.c, $(SRCS_NOMAIN))
OBJS_NOMAIN := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS_NOMAIN))
OBJS := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))


.PHONY: all clean

# library and several programs
TARGET = libgrt grt grt.syn grt.travt

all: $(TARGET) 


libgrt:      $(LIB_DIR)/$(LIB_NAME)
grt:         $(BIN_DIR)/grt 
grt.syn:     $(BIN_DIR)/grt.syn 
grt.travt:	 $(BIN_DIR)/grt.travt


# generate executable files
$(BIN_DIR)/grt: $(BUILD_DIR)/grt_main.o  libgrt
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(CFLAGS) $(FOMPFLAGS) $(RPATH_FLAGS)

$(BIN_DIR)/grt.syn: $(BUILD_DIR)/grt_syn.o  libgrt
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(CFLAGS) $(FOMPFLAGS) $(RPATH_FLAGS)

$(BIN_DIR)/grt.travt: $(BUILD_DIR)/grt_travt.o  libgrt
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(CFLAGS) $(FOMPFLAGS) $(RPATH_FLAGS)


# Link dynamic library
# The order in the compiled statement is critical, otherwise the library will not be linked
$(LIB_DIR)/$(LIB_NAME): $(OBJS_NOMAIN)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(LIB_DIR)
	$(CC) -o $@ $^ $(LDFLAGS)  $(FOMPFLAGS) 
	

# Compile object file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC)  -o $@ -c $< $(CFLAGS) $(FOMPFLAGS)


# Automatically generate dependencies
DEPENDS := $(OBJS:.o=.d)
-include $(DEPENDS)

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	@$(CC) $(CFLAGS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,$(BUILD_DIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# remove all generated files
clean:
	rm -rf $(BIN_DIR)
	rm -rf $(BUILD_DIR)
	rm -rf $(LIB_DIR) 