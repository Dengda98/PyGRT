BUILD_DIR := build
INC_DIR := include
SRC_DIR := src
BIN_DIR := bin
LIB_DIR := lib
LNAME := grt
LIBS_NAME := lib$(LNAME).a
LIBD_NAME := lib$(LNAME).so

VERSION := v$(shell cat version )


# link options for FFTW
# NEED EXPLICITY SET LIBRARY PATH on MacOS
LFFTW_FLAGS := -l:libfftw3.a  -l:libfftw3f.a

CC := gcc


# add -static to enforce gcc link the static library, no matter standard-lib or custom-lib,
# -l like "-l:libXYZ.a" should be saved for conditional compile.
FOMPFLAGS := -fopenmp 

# NetCDF library
LNETCDF_FLAGS :=  -lnetcdf

# link static library on Windows
LINK_STATIC := 

# expand stack memory for Windows
STACK_MEM := 

ifeq ($(OS),Windows_NT)  # link static oenpmp on windows
	STACK_MEM := -Wl,-stack,0x1000000
    LINK_STATIC := -static
endif

LDFLAGS := $(LINK_STATIC) $(FOMPFLAGS) -lm $(LFFTW_FLAGS) $(LNETCDF_FLAGS) $(STACK_MEM)

# change architecture for macOS, from make command
ARCH = 

# load args from command
CFLAGS2 := 

# However, Maybe -static is not working for standard-lib on MacOS, 
# at least for now, it's not a big problem.
CFLAGS := -O3 -std=gnu99 \
		  -fPIC -Wall -Wextra -I$(shell realpath $(INC_DIR)) $(CFLAGS2) \
		  -DGRT_VERSION=\"$(VERSION)\" -D_GNU_SOURCE  $(ARCH)  $(FOMPFLAGS)
#  -fdump-tree-all  -g -ffast-math -O3 -fno-associative-math  -march=native -mtune=native 
#  -Wconversion -Wsign-conversion -Wfloat-conversion -Wcast-qual -Wcast-align


# 获取所有源文件
SRCS = $(shell find $(SRC_DIR) -name '*.c')
# 转为中间文件
OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))

DEPS := $(OBJS:.o=.d)  # include main functions here

PROGS := $(BIN_DIR)/grt

# Main targets
.PHONY: all objs progs libs clean cleanbuild

all: objs progs libs

objs: $(BUILD_DIR) $(OBJS)

progs: objs $(BIN_DIR) $(PROGS)

libs: objs $(LIB_DIR) $(LIB_DIR)/$(LIBS_NAME)  $(LIB_DIR)/$(LIBD_NAME)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) -o $@  -c $< $(CFLAGS) 

# ------------------------Executable files------------------------
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

$(PROGS): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# ----------------------- Dependency generation -----------------------
# Dependency generation (only for non-clean targets)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleanbuild)
-include $(DEPS)
endif
endif

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.c
	@mkdir -p $(shell dirname $@)
	@$(CC) $(CFLAGS) -MM $< > $@.$$$$; \
	sed 's,\($(notdir $*)\)\.o[ :]*,$(shell dirname $@)/\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# ----------------------------------------------------------------------

$(LIB_DIR):
	@mkdir -p $@

# build dynamic library
# The order in the compiled statement is critical, otherwise the library will not be linked
$(LIB_DIR)/$(LIBD_NAME): $(OBJS)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

# build static library
$(LIB_DIR)/$(LIBS_NAME): $(OBJS)
	ar rcs $@ $^  


cleanbuild:
	rm -rf $(BUILD_DIR)

clean: cleanbuild
	rm -rf $(LIB_DIR)
	rm -rf $(BIN_DIR)