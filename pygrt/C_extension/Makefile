BUILD_DIR := build
INC_DIR := include
BIN_DIR := bin
LIB_DIR := lib
LNAME := grt
LIBS_NAME := lib$(LNAME).a
LIBD_NAME := lib$(LNAME).so

VERSION := v$(shell cat version )


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
	STACK_MEM := -Wl,-stack,0x1000000
    LINK_STATIC := -static
	LDFLAGS := $(LINK_STATIC) $(FOMPFLAGS)
endif

# change architecture for macOS, from make command
ARCH = 

# However, Maybe -static is not working for standard-lib on MacOS, 
# at least for now, it's not a big problem.
CFLAGS :=  $(LINK_STATIC) $(LFFT_FLAGS) -lm -O3 \
		-fPIC -Wall -Wextra $(STACK_MEM) -I$(shell realpath $(INC_DIR)) \
		  -D'GRT_VERSION=\"$(VERSION)\"'  $(ARCH)  $(FOMPFLAGS)
#  -fdump-tree-all  -g -ffast-math -O3 -fno-associative-math  -march=native -mtune=native 


# 定义子目录
SUBDIRS = \
    src/common \
	src/dynamic \
	src/static \
	src/travt


all: objs progs libs


objs:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir CC="$(CC)" CFLAGS="$(CFLAGS)" objs || exit 1; \
	done

progs: objs 
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir CC="$(CC)" CFLAGS="$(CFLAGS)" progs || exit 1; \
	done


libs: objs $(LIB_DIR) $(LIB_DIR)/$(LIBS_NAME)  $(LIB_DIR)/$(LIBD_NAME)

$(LIB_DIR):
	@mkdir -p $@

# build dynamic library
# The order in the compiled statement is critical, otherwise the library will not be linked
OBJS = $(wildcard $(BUILD_DIR)/**/*.o)
$(LIB_DIR)/$(LIBD_NAME): objs
	$(CC) -shared -o $@ $(OBJS) $(LDFLAGS) 

# build static library
$(LIB_DIR)/$(LIBS_NAME): objs
	ar rcs $@ $(OBJS)  
	

cleanbuild:
	rm -rf $(BUILD_DIR)

clean: cleanbuild
	rm -rf $(LIB_DIR)
	rm -rf $(BIN_DIR)