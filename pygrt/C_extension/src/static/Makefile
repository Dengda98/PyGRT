
CC := gcc 
CFLAGS := -Wall -g -fPIC -I../../include -lm

COMMON_OBJS := $(wildcard ../../build/common/*.o)

BUILD_DIR = ../../build/dynamic
BIN_DIR = ../../bin
SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRCS))

DEPS := $(OBJS:.o=.d)  # include main functions here

OBJS := $(filter-out \
	$(BUILD_DIR)/stgrt_main.o \
	$(BUILD_DIR)/stgrt_syn.o \
	$(BUILD_DIR)/stgrt_strain.o \
	$(BUILD_DIR)/stgrt_stress.o \
	$(BUILD_DIR)/stgrt_rotation.o \
	, $(OBJS))

PROGS := $(BIN_DIR)/stgrt \
		 $(BIN_DIR)/stgrt.syn \
		 $(BIN_DIR)/stgrt.strain \
		 $(BIN_DIR)/stgrt.stress \
		 $(BIN_DIR)/stgrt.rotation

all: objs progs

objs: $(BUILD_DIR) $(OBJS)

progs: $(OBJS) $(BIN_DIR) $(PROGS)


$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.c
	$(CC) -o $@  -c $< $(CFLAGS) 


# ------------------------Executable files------------------------
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

$(BIN_DIR)/stgrt: stgrt_main.c $(OBJS) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

$(BIN_DIR)/stgrt.syn: stgrt_syn.c $(OBJS) $(COMMON_OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

$(BIN_DIR)/stgrt.strain: stgrt_strain.c $(COMMON_OBJS) 
	$(CC) -o $@ $^ $(CFLAGS)

$(BIN_DIR)/stgrt.stress: stgrt_stress.c $(COMMON_OBJS) 
	$(CC) -o $@ $^ $(CFLAGS)

$(BIN_DIR)/stgrt.rotation: stgrt_rotation.c $(COMMON_OBJS) 
	$(CC) -o $@ $^ $(CFLAGS)

# ----------------------- Dependency generation -----------------------
-include $(DEPS)

$(BUILD_DIR)/%.d: %.c
	@mkdir -p $(shell dirname $@)
	@$(CC) $(CFLAGS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,$(BUILD_DIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# ---------------------------------------------------------------------

cleanbuild:
	rm -rf $(BUILD_DIR)

clean: cleanbuild
	rm -f $(PROGS)