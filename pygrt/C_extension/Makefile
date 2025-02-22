
SRC_DIR := src
INC_DIR := include
BUILD_DIR := build
BIN_DIR := bin
LIB_DIR := lib
LNAME := grt
LIBS_NAME := lib$(LNAME).a
LIBD_NAME := lib$(LNAME).so

VERSION := v$(shell cat version )
CURRENT_DIR := $(shell pwd)

# use static library, make life easier...
# RPATH_FLAGS :=  -L$(LIB_DIR) -l:$(LIBS_NAME)  # ':libXYZ.a' not worked on MacOS
RPATH_FLAGS := $(LIB_DIR)/$(LIBS_NAME)

# link options for FFTW
# NEED EXPLICITY SET LIBRARY PATH on MacOS
LFFT_FLAGS := -l:libfftw3.a  -l:libfftw3f.a

CC := gcc
FOMPFLAGS := -fopenmp 

CFLAGS := $(LFFT_FLAGS) -O3 -g \
          -fPIC -Wall -Wextra -I$(INC_DIR) \
		  -lm -DGRT_VERSION="\"$(VERSION)\""
#  -fdump-tree-all  -g -ffast-math -O3 -fno-associative-math  -march=native -mtune=native 

# change architecture for macOS
ifdef ARCH
  CFLAGS += -arch $(ARCH)
endif

# options for dynamic library
LDFLAGS := -shared -lm


INCS := $(wildcard $(INC_DIR)/*.h)
SRCS := $(wildcard $(SRC_DIR)/*.c)
SRCS_NOMAIN := $(SRCS)

# remove main functions
SRCS_NOMAIN := $(filter-out $(SRC_DIR)/grt_main.c, $(SRCS_NOMAIN))
SRCS_NOMAIN := $(filter-out $(SRC_DIR)/grt_syn.c, $(SRCS_NOMAIN))
SRCS_NOMAIN := $(filter-out $(SRC_DIR)/grt_travt.c, $(SRCS_NOMAIN))
OBJS_NOMAIN := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS_NOMAIN))
OBJS := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))


.PHONY: all clean cleanbuild

# library and several programs
TARGET = libgrt grt grt.syn  grt.travt

all: $(TARGET) 


libgrt:      $(LIB_DIR)/$(LIBS_NAME)  $(LIB_DIR)/$(LIBD_NAME)
grt:         $(BIN_DIR)/grt 
grt.syn:     $(BIN_DIR)/grt.syn 
grt.travt:	 $(BIN_DIR)/grt.travt


# generate executable files
$(BIN_DIR)/grt: $(BUILD_DIR)/grt_main.o  libgrt
	@echo $(shell ls $(LIB_DIR)/*)
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(RPATH_FLAGS) $(CFLAGS) $(FOMPFLAGS)

$(BIN_DIR)/grt.syn: $(BUILD_DIR)/grt_syn.o  libgrt
	@echo $(shell ls $(LIB_DIR)/*)
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(RPATH_FLAGS) $(CFLAGS) $(FOMPFLAGS)

$(BIN_DIR)/grt.travt: $(BUILD_DIR)/grt_travt.o  libgrt
	@echo $(shell ls $(LIB_DIR)/*)
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(RPATH_FLAGS) $(CFLAGS) $(FOMPFLAGS)


# build dynamic library
# The order in the compiled statement is critical, otherwise the library will not be linked
$(LIB_DIR)/$(LIBD_NAME): $(OBJS_NOMAIN)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(LIB_DIR)
	$(CC) -o $@ $^ $(LDFLAGS)  $(FOMPFLAGS) 

# build static library
$(LIB_DIR)/$(LIBS_NAME): $(OBJS_NOMAIN)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(LIB_DIR)
	ar rcs $@ $^ 
	

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

# remove object files
cleanbuild:
	rm -rf $(BUILD_DIR)