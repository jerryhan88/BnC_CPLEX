SYSTEM:=x86-64_osx
DYL_EXT:=dylib

CXX:=clang++
C:=clang

CPLEX_LN_FLAGS:=-lconcert -lilocplex -lcplex
CPLEX_LN_FLAGS+=-m64 -lm -lpthread 
CPLEX_LN_FLAGS+=-framework CoreFoundation 
CPLEX_LN_FLAGS+=-framework IOKit
CPLEX_LN_FLAGS+=-stdlib=libc++

CPLEX_HOME:=/Applications/CPLEX_Studio128
UTIL_HOME:=/Users/ckhan/workspace/util
JSON_HOME:=/Users/ckhan/workspace/json-develop

SRC_DIR:=./src
BUILD_DIR:=./build
LIB_DIR:=./lib

LIB_FORMAT:=static_pic
SRC_EXT:=cpp
OBJ_EXT:=o

CXX_OPT:=-std=c++11 -DIL_STD
#####################################################################
#
# Set paths
#
CPLEX_DIR:=$(CPLEX_HOME)/cplex
CPLEX_INC_DIR:=$(CPLEX_DIR)/include
CPLEX_LIB_DIR:=$(CPLEX_DIR)/lib/$(SYSTEM)/$(LIB_FORMAT)

CONCERT_DIR:=$(CPLEX_HOME)/concert
CONCERT_INC_DIR:=$(CONCERT_DIR)/include
CONCERT_LIB_DIR:=$(CONCERT_DIR)/lib/$(SYSTEM)/$(LIB_FORMAT)

UTIL_NAME:=util
UTIL_INC_DIR:=$(UTIL_HOME)/include
UTIL_LIB_DIR:=$(UTIL_HOME)/lib

JSON_INC_DIR:=$(JSON_HOME)/single_include

#####################################################################
MAIN_CPP:=main.$(SRC_EXT)
MAIN_OBJ:=$(BUILD_DIR)/main.$(OBJ_EXT)

#####################################################################
SOURCES:=$(shell find $(SRC_DIR) -type f -name "*.$(SRC_EXT)")
OBJECTS:=$(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.$(SRC_EXT)=.$(OBJ_EXT)))

#####################################################################
SEQ_NAME:=sequencer
SEQ_DIR:=$(BUILD_DIR)/GH
SEQ_C:=$(SRC_DIR)/GH/Sequencer.c
SEQ_OBJ:=$(SEQ_DIR)/Sequencer.$(OBJ_EXT)
SEQ_STL:=$(LIB_DIR)/lib$(SEQ_NAME).a
SEQ_DYL:=$(LIB_DIR)/lib$(SEQ_NAME).$(DYL_EXT)

#####################################################################
CXX_INC_DIRS:=-I$(CPLEX_INC_DIR) -I$(CONCERT_INC_DIR) -I$(UTIL_INC_DIR)
CXX_LIB_DIRS:=-L$(CPLEX_LIB_DIR) -L$(CONCERT_LIB_DIR) -L$(UTIL_LIB_DIR) -L$(LIB_DIR)

CXX_FLAGS:=$(CXX_OPT) $(CXX_INC_DIRS)
CXX_LN_FLAGS:=$(CPLEX_LN_FLAGS) -l$(SEQ_NAME) -l$(UTIL_NAME)

#####################################################################
PROGRAM:=router
PROGRAM_STL:=$(LIB_DIR)/lib$(PROGRAM).a

#####################################################################

all: $(PROGRAM) $(PROGRAM_STL)

echoTest:
	@echo "echo TEST"
	@echo $(SOURCES)

#Target:Pre-req
$(PROGRAM): $(MAIN_OBJ) $(OBJECTS) $(SEQ_STL)
	$(CXX) $(CXX_FLAGS) $(CXX_LIB_DIRS) $(MAIN_OBJ) $(OBJECTS) -o $(PROGRAM) $(CXX_LN_FLAGS) 

$(MAIN_OBJ): $(MAIN_CPP)
	@mkdir -p $(dir $@)
	@echo "Start compile"
	$(CXX) -c $(CXX_FLAGS) $< -o $@
	@echo "Build Success!!!"
	@echo ""

#Compile each .cpp file
$(BUILD_DIR)/%.$(OBJ_EXT): $(SRC_DIR)/%.$(SRC_EXT)
	@mkdir -p $(dir $@)
	@echo "Start compile"
	@echo "Target file:" $@ "; Pre-req:" $<
	$(CXX) -c $(CXX_FLAGS) $< -o $@
	@echo "Build Success!!!"
	@echo ""

$(PROGRAM_STL): $(OBJECTS)
	@mkdir -p $(dir $@)
	@echo "Build a router static library"	
	ar rvs $(PROGRAM_STL) $(OBJECTS)
	@echo "Successfully build the library"	
	@echo ""

$(SEQ_STL): $(SEQ_C)
	@mkdir -p $(dir $@)
	@echo "Build a sequencer static library"	
	$(C) -c $(SEQ_C) -o $(SEQ_OBJ)
	ar rvs $(SEQ_STL) $(SEQ_OBJ)
	@echo "Successfully build the library"	
	@echo ""

$(SEQ_DYL): $(SEQ_C)
	@mkdir -p $(dir $@)
	@echo "Build a sequencer dynamic library"	
	$(C) -c -fPIC $(SEQ_C) -o $(SEQ_OBJ)
	$(C) -shared $(SEQ_OBJ) -o $(SEQ_DYL)
	@echo "Successfully build the library"	
	@echo ""

clean:
	@rm -rf $(BUILD_DIR)
	@rm -rf $(LIB_DIR)
	@rm $(PROGRAM)