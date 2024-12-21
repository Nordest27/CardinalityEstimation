# Compiler
CXX = g++

# Compiler Flags
CXXFLAGS = -std=c++17 -Wall -O3 -march=native

# Target executable
TARGET = cardinality_estimation.bin

# Source files
SRCS = main.cc

# Header files
HDRS = xxHash64.h

# Build target
$(TARGET): $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

# Clean target
clean:
	rm -f $(TARGET)
