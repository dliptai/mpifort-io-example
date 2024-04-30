# Compiler
FC = mpifort
DEBUG = no

# Compiler flags
FCFLAGS =

ifeq ($(DEBUG),yes)
FCFLAGS += -Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
endif

# Target executable
TARGET = main

# Source files
SRCS = utils.f90 io_write.f90 io_read.f90 main.f90

# Object files
OBJS = $(SRCS:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f $(OBJS) $(TARGET) $(SRCS:.f90=.mod)
