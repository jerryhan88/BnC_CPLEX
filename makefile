SYSTEM := x86-64_osx
LIB_FORMAT := static_pic

CPLEX_HOME := /Applications/CPLEX_Studio128

SRC_DIR := /Users/ckhan/workspace/BnC_CPLEX/BnC_CPLEX
BUILD_DIR := ./build

CXX_LN_FLAGS := -lconcert -lilocplex -lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit -stdlib=libc++

#####################################################################

CPLEX_DIR := $(CPLEX_HOME)/cplex
CPLEX_INC_DIR := $(CPLEX_DIR)/include
CPLEX_LIB_DIR := $(CPLEX_DIR)/lib/$(SYSTEM)/$(LIB_FORMAT)

CONCERT_DIR := $(CPLEX_HOME)/concert
CONCERT_INC_DIR := $(CONCERT_DIR)/include
CONCERT_LIB_DIR := $(CONCERT_DIR)/lib/$(SYSTEM)/$(LIB_FORMAT)

CXX = clang++
CXX_OPT := -std=c++11 -DIL_STD
CXX_FLAGS := $(CXX_OPT) -I$(CPLEX_INC_DIR) -I$(CONCERT_INC_DIR)
CXX_LN_DIRS := -L$(CPLEX_LIB_DIR) -L$(CONCERT_LIB_DIR)

SRC_EXT := cpp
OBJ_EXT := o
SOURCES := $(shell find $(SRC_DIR) -type f -name *.$(SRC_EXT))
OBJECTS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.$(SRC_EXT)=.$(OBJ_EXT)))


C = clang
DYL_EXT := dylib

SEQ_NAME := Sequencer

SEQ_DIR := $(BUILD_DIR)/GH
SEQ_C := $(SRC_DIR)/GH/$(SEQ_NAME).c
SEQ_O := $(SEQ_DIR)/$(SEQ_NAME).o
SEQ_D := $(SEQ_DIR)/lib$(SEQ_NAME).$(DYL_EXT)

CXX_LN_FLAGS += -l$(SEQ_NAME)
CXX_LN_DIRS += -L$(SEQ_DIR)

PROGRAM := bnc

all: $(PROGRAM)

echoTest:
	@echo "echo TEST"
	@echo $(OBJECTS)

#Target:Pre-req
$(PROGRAM): $(OBJECTS) $(SEQ_D)
	$(CXX) $(CXX_FLAGS) $(CXX_LN_DIRS) $(OBJECTS) -o $(PROGRAM) $(CXX_LN_FLAGS)

#Compile each .cpp file
$(BUILD_DIR)/%.$(OBJ_EXT): $(SRC_DIR)/%.$(SRC_EXT)
	@mkdir -p $(dir $@)
	@echo "Compile a object file"
	@echo "Target file:" $@ "; Pre-req:" $<
	$(CXX) -c $(CXX_FLAGS) $< -o $@
	@echo "Build Success!!!"
	@echo ""

$(SEQ_D): $(SEQ_C)
	@mkdir -p $(dir $@)
	@echo "Build sequencer"	
	$(C) -c -fPIC $(SEQ_C) -o $(SEQ_O)
	$(C) -shared $(SEQ_O) -o $(SEQ_D)
	@echo ""

clean:
	@rm -rf $(BUILD_DIR)
	@rm $(PROGRAM)