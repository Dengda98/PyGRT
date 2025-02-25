
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


# add -static to enforce gcc link the static library, no matter standard-lib or custom-lib,
# -l like "-l:libXYZ.a" should be saved for conditional compile.
FOMPFLAGS := -fopenmp 

# link static library on Windows
LINK_STATIC := 
LDFLAGS := $(FOMPFLAGS)

# expand stack memory for Windows
STACK_MEM := 

ifeq ($(OS),Windows_NT)  # link static oenpmp on windows
	STACK_MEM := -Wl,-stack_size,0x1000000
    LINK_STATIC := -static
	LDFLAGS := $(LINK_STATIC) $(FOMPFLAGS)
endif

# change architecture for macOS, from make command
ARCH = 

# However, Maybe -static is not working for standard-lib on MacOS, 
# at least for now, it's not a big problem.
CFLAGS :=  $(LINK_STATIC) $(LFFT_FLAGS) -lm -O3 -g \
		-fPIC -Wall -Wextra $(STACK_MEM) -I$(INC_DIR) \
		  -DGRT_VERSION="\"$(VERSION)\""  $(ARCH)  $(FOMPFLAGS)
#  -fdump-tree-all  -g -ffast-math -O3 -fno-associative-math  -march=native -mtune=native 


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
TARGET = libgrt \
         grt \
		 grt.syn \
		 grt.travt

all: $(TARGET) 


libgrt:      $(LIB_DIR)/$(LIBS_NAME)  $(LIB_DIR)/$(LIBD_NAME)
grt:         $(BIN_DIR)/grt 
grt.syn:     $(BIN_DIR)/grt.syn 
grt.travt:	 $(BIN_DIR)/grt.travt


# generate executable files
$(BIN_DIR)/grt: $(BUILD_DIR)/grt_main.o  libgrt
	@echo $(shell ls $(LIB_DIR)/*)
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(RPATH_FLAGS) $(CFLAGS)

$(BIN_DIR)/grt.syn: $(BUILD_DIR)/grt_syn.o  libgrt
	@echo $(shell ls $(LIB_DIR)/*)
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(RPATH_FLAGS) $(CFLAGS)

$(BIN_DIR)/grt.travt: $(BUILD_DIR)/grt_travt.o  libgrt
	@echo $(shell ls $(LIB_DIR)/*)
	@mkdir -p $(BIN_DIR)
	$(CC) -o $@ $<  $(RPATH_FLAGS) $(CFLAGS)


# build dynamic library
# The order in the compiled statement is critical, otherwise the library will not be linked
$(LIB_DIR)/$(LIBD_NAME): $(OBJS_NOMAIN)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(LIB_DIR)
	$(CC) -shared -o $@ $^ $(LDFLAGS) 

# build static library
$(LIB_DIR)/$(LIBS_NAME): $(OBJS_NOMAIN)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(LIB_DIR)
	ar rcs $@ $^ 
	

# Compile object file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC)  -o $@ -c $< $(CFLAGS)


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