
CC := gcc 
CFLAGS := -Wall -g -fPIC -I../../include -lm

BUILD_DIR = ../../build/common
SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

progs: 
	@echo "No rule to make target 'progs'"

all: objs

objs: $(BUILD_DIR) $(OBJS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)


$(BUILD_DIR)/%.o: %.c
	$(CC) -o $@  -c $< $(CFLAGS) 


# ----------------------- Dependency generation -----------------------
-include $(DEPS)

$(BUILD_DIR)/%.d: %.c
	@mkdir -p $(shell dirname $@)
	@$(CC) $(CFLAGS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,$(BUILD_DIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# ---------------------------------------------------------------------

cleanbuild: clean

clean:
	rm -rf $(BUILD_DIR)