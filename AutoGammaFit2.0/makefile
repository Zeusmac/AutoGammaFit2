
INCLUDES := \
    -I. \
	-I./app \
    -I./core \
    -I./fitting \
    -I./tracking \
    -I./models \
    -I./database \
    -I./io \
    -I./analysis_modules \
    -I./gui

# =========================
# ROOT CONFIG
# =========================
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs) -lSpectrum

CXX         := g++
CXXFLAGS = -O2 -std=c++17 -Wall -fPIC -pthread -m64 \
           $(ROOTCFLAGS) \
           $(INCLUDES)
LDFLAGS     := $(ROOTLIBS)

# =========================
# DIRECTORIES
# =========================
SRC_DIRS := app core fitting tracking models database io analysis_modules
OBJ_DIR  := build
BIN_DIR   := bin

TARGET   := $(BIN_DIR)/gamma_fit

# =========================
# SOURCES
# =========================
SOURCES := $(shell find $(SRC_DIRS) -name "*.cpp")
OBJECTS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# =========================
# DEFAULT RULE
# =========================
all: $(TARGET)

# =========================
# LINK
# =========================
$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# =========================
# COMPILE
# =========================
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# =========================
# CREATE DIRS
# =========================
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# =========================
# CLEAN
# =========================
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# =========================
# RUN
# =========================
run: all
	./$(TARGET)

# =========================
# GUI TARGET
# =========================
GUI_DIR      := gui
GUI_TARGET   := $(BIN_DIR)/gamma_gui
GUI_ROOTLIBS := $(shell root-config --libs) -lSpectrum -lGui -lGuiHtml

# Core objects: everything except the CLI main
CORE_SRCS := $(filter-out app/FitGammaAnalyzer.cpp, \
                 $(shell find core fitting tracking models database io analysis_modules -name "*.cpp"))
CORE_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(CORE_SRCS))

# GUI-specific sources
GUI_SRCS     := $(GUI_DIR)/GammaFitGUI.cpp $(GUI_DIR)/RunGUI.cpp
GUI_OWN_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(GUI_SRCS))

# rootcling dictionary
GUI_DICT_SRC := $(OBJ_DIR)/gui/GammaFitGUIDict.cpp
GUI_DICT_OBJ := $(OBJ_DIR)/gui/GammaFitGUIDict.o

$(GUI_DICT_SRC): $(GUI_DIR)/GammaFitGUI.h $(GUI_DIR)/LinkDef.h
	@mkdir -p $(dir $@)
	rootcling -f $@ -c $(INCLUDES) -I$(shell root-config --incdir) $^

$(GUI_DICT_OBJ): $(GUI_DICT_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/gui/%.o: $(GUI_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: gui
gui: $(GUI_TARGET)

$(GUI_TARGET): $(CORE_OBJS) $(GUI_OWN_OBJS) $(GUI_DICT_OBJ) | $(BIN_DIR)
	$(CXX) $^ -o $@ $(GUI_ROOTLIBS)