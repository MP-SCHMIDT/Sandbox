TARGET_EXEC = main.exe

BUILD_DIR = ./Build

SRCS = $(wildcard Source/*.cpp) $(wildcard Source/Algo/*/*.cpp) $(wildcard Source/Projects/*/*.cpp)
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
DEPS = $(SRCS:%.cpp=$(BUILD_DIR)/%.d)

CXX = g++ -std=c++23

FLAGS_DEPEND = -MMD -MP
FLAGS_ENV = -m64
FLAGS_PARALLEL = -fopenmp -pthread -lOpenCL
FLAGS_WARNING = -Wall -Wextra
# FLAGS_DEBUG = -g
# FLAGS_OPTIM =
FLAGS_DEBUG = -g1
FLAGS_OPTIM = -O2 -march=native
CXXFLAGS = $(FLAGS_DEPEND) $(FLAGS_ENV) $(FLAGS_PARALLEL) $(FLAGS_WARNING) $(FLAGS_DEBUG) $(FLAGS_OPTIM)

UNAME := $(shell uname)
ifeq ($(UNAME), Windows_NT)
LIB_FLAGS_GL = -lfreeglut -lopengl32 -lglu32
else ifeq ($(UNAME), MINGW32_NT-6.2)
LIB_FLAGS_GL = -lfreeglut -lopengl32 -lglu32
else
LIB_FLAGS_GL = -lGL -lglut -lGLU -lX11 -lm
endif

LIB_FLAGS = -L"Libs/freeglut/lib/x64" -L"Libs/OpenCL/lib" $(LIB_FLAGS_GL)
INC_FLAGS = -I"Libs" -I"Libs/freeglut/include" -I"Libs/OpenCL/include" -I"Source" -I"Source/Algo" -I"Source/Projects"

$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LIB_FLAGS) $(CXXFLAGS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(INC_FLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET_EXEC)

-include $(DEPS)
