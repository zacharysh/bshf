TARGET_EXEC := bshf

SRC_DIR  = ./src
BUILD_DIR= ./build

SRCS := $(shell find $(SRC_DIR) -type f -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

INC_DIRS := $(shell find $(SRC_DIR) -type d) 
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CC=g++ -march=native -fopenmp
CPPFLAGS= -std=c++17 -Werror -pedantic -Wextra -Wdouble-promotion -Wconversion -O3
LDFLAGS= -llapack -lblas 

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
