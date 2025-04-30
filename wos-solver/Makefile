# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wpedantic -Wunused -Wunreachable-code -Werror

# Source files
SRCS = main.cpp wost_algorithm.cpp experiment_runner.cpp csv_io.cpp

# Object files (optional for larger projects, here we stay simple)
OBJS = $(SRCS:.cpp=.o)

# Output binary name
TARGET = wos

# Default build target
all: $(TARGET)

# How to build the target
$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

# Clean build files
clean:
	rm -f $(TARGET) $(OBJS)
