
CC := gcc 
CFLAGS := -Wall -g -fPIC -I../../include -lm

COMMON_OBJS := $(wildcard ../../build/common/*.o)
TRAVT_OBJS := $(wildcard ../../build/travt/*.o)

BUILD_DIR = ../../build/dynamic
BIN_DIR = ../../bin
SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRCS))

DEPS := $(OBJS:.o=.d)  # include main functions here

OBJS := $(filter-out \
	$(BUILD_DIR)/grt_main.o \
	$(BUILD_DIR)/grt_syn.o \
	$(BUILD_DIR)/grt_b2a.o \
	$(BUILD_DIR)/grt_strain.o \
	$(BUILD_DIR)/grt_stress.o \
	$(BUILD_DIR)/grt_rotation.o \
	, $(OBJS))

PROGS := $(BIN_DIR)/grt \
         $(BIN_DIR)/grt.syn \
		 $(BIN_DIR)/grt.b2a \
		 $(BIN_DIR)/grt.strain \
		 $(BIN_DIR)/grt.stress \
		 $(BIN_DIR)/grt.rotation

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

$(BIN_DIR)/grt: grt_main.c $(OBJS) $(COMMON_OBJS) $(TRAVT_OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

$(BIN_DIR)/grt.syn: grt_syn.c $(OBJS) $(COMMON_OBJS) $(TRAVT_OBJS)
	$(CC) -o $@ $^  $(CFLAGS)

$(BIN_DIR)/grt.b2a: grt_b2a.c $(COMMON_OBJS)
	$(CC) -o $@ $^  $(CFLAGS)

$(BIN_DIR)/grt.strain: grt_strain.c $(COMMON_OBJS)
	$(CC) -o $@ $^  $(CFLAGS)

$(BIN_DIR)/grt.stress: grt_stress.c $(COMMON_OBJS)
	$(CC) -o $@ $^  $(CFLAGS)

$(BIN_DIR)/grt.rotation: grt_rotation.c $(COMMON_OBJS)
	$(CC) -o $@ $^  $(CFLAGS)

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