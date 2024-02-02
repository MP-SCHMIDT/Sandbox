TARGET_EXEC = main.exe

BUILD_DIR = ./Build
SRC_DIRS = ./Source

SRCS = $(wildcard Source/*.cpp) $(wildcard Source/*/*.cpp) $(wildcard Source/*/*/*.cpp)

OBJS = $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS = $(SRCS:%=$(BUILD_DIR)/%.d)

CXX = g++ -std=c++23

FLAGS_ENV = -m64
FLAGS_OPTIM = -O2 -march=native
FLAGS_OPENMP = -fopenmp
FLAGS_DEBUG = -g1
FLAGS_WARNING = -Wall -Wextra
FLAGS_DEPEND = -MMD -MP
CXXFLAGS = $(FLAGS_DEPEND) $(FLAGS_ENV) $(FLAGS_OPENMP) $(FLAGS_WARNING) $(FLAGS_DEBUG) $(FLAGS_OPTIM)

UNAME := $(shell uname)
ifeq ($(UNAME), Windows_NT)
FLAGS_GL = -lfreeglut -lopengl32 -lglu32
else ifeq ($(UNAME), MINGW32_NT-6.2)
FLAGS_GL = -lfreeglut -lopengl32 -lglu32
else
FLAGS_GL = -lGL -lglut -lGLU -lX11 -lm
endif

LIB_FLAGS = -L"Libs/freeglut/lib/x64" $(FLAGS_GL)
INC_FLAGS = -I"Libs" -I"Source" -I"Source/Algo" -I"Source/Projects"

$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LIB_FLAGS) $(CXXFLAGS)

$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(INC_FLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET_EXEC)

-include $(DEPS)
