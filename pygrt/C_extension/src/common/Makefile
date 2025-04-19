
CC := gcc 
CFLAGS := -Wall -g -fPIC -I../../include -lm

BUILD_DIR = ../../build/common
BIN_DIR = ../../bin
SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

OBJS := $(filter-out \
	$(BUILD_DIR)/grt_k2a.o \
	, $(OBJS))

PROGS := $(BIN_DIR)/grt.k2a

progs: $(OBJS) $(BIN_DIR) $(PROGS)

all: objs progs

objs: $(BUILD_DIR) $(OBJS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)


$(BUILD_DIR)/%.o: %.c
	$(CC) -o $@  -c $< $(CFLAGS) 

# ------------------------Executable files------------------------
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

$(BIN_DIR)/grt.k2a: grt_k2a.c $(OBJS)
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