TARGET_EXEC := schrodingerHF

SRC_DIR  = ./src
BUILD_DIR= ./build

SRCS := $(shell find $(SRC_DIR) -type f -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

INC_DIRS := $(shell find $(SRC_DIR) -type d) 
INC_FLAGS := $(addprefix -I,$(INC_DIRS)) -I/opt/homebrew/opt/lapack/include

LDFLAGS=-L/opt/homebrew/opt/lapack/lib -llapack -lblas

CC=g++
CPPFLAGS= -O2 -fno-inline -Werror -pedantic -Wextra -Wdouble-promotion -Wconversion
# -Wshadow -Wsign-conversion


# The final build step.
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)
	mv $@ .

# Build step for C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean:
	@echo "Cleaning up..."
	@rm -rf $(BUILD_DIR)
