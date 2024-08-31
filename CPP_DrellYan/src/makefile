# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -std=c++17 -I.

# Source files
SRCS = main.cpp functions.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Output executable name
TARGET = my_program


# Linker flags
LDFLAGS = -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2/build/bin -lcuba

# Build the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean the build
clean:
	rm -f $(OBJS) $(TARGET)

