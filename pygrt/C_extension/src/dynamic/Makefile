
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
	, $(OBJS))

PROGS := $(BIN_DIR)/grt \
         $(BIN_DIR)/grt.syn

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

$(BIN_DIR)/grt: grt_main.c $(OBJS)
	$(CC) -o $@ $^  $(COMMON_OBJS) $(TRAVT_OBJS)  $(CFLAGS)

$(BIN_DIR)/grt.syn: grt_syn.c $(OBJS)
	$(CC) -o $@ $^  $(COMMON_OBJS) $(TRAVT_OBJS)  $(CFLAGS)

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