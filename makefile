# Define the compiler
CXX = g++
# Define the target executable
TARGET = main
# Define source files
SRCS = main.cpp func.cpp io.cpp common.cpp export.cpp
# Object directory
OBJDIR = ../build
# Define object files
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)

# Default rule
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Compile each .cpp file into .o files in OBJDIR
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure the build directory exists
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Clean rule to remove object files and the executable with the folder
clean:
	rm -rf $(OBJDIR) $(TARGET)

# Phony targets
.PHONY: all clean