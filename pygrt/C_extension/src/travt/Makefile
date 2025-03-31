
CC := gcc 
CFLAGS := -Wall -g -fPIC -I../../include  -lm

COMMON_OBJS := $(wildcard ../../build/common/*.o)

BUILD_DIR = ../../build/travt
BIN_DIR = ../../bin
SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRCS))

DEPS := $(OBJS:.o=.d)  # include main functions here

OBJS := $(filter-out \
	$(BUILD_DIR)/grt_travt.o \
	, $(OBJS))

PROGS := $(BIN_DIR)/grt.travt

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

$(BIN_DIR)/grt.travt: grt_travt.c $(OBJS) $(COMMON_OBJS) 
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