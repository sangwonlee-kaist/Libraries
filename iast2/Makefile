CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall
SRCS = $(wildcard *.cpp)
OBJS_ALL = $(SRCS:.cpp=.o)
OBJS = $(patsubst arma_solver.%, , $(OBJS_ALL)) # remove arma solver. we don't need it.
TARGET = libiast.a

#.SUFFIXES : .cpp. o

all: $(TARGET)

$(TARGET): $(OBJS)
	ar rc $(@) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $(<)

clean:
	$(RM) $(OBJS) $(TARGET)
