
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
# GUI TARGET NAMES
# =========================
GUI_DIR      := gui
GUI_TARGET   := $(BIN_DIR)/gamma_gui
GUI_ROOTLIBS := $(shell root-config --libs) -lSpectrum -lGui -lGuiHtml -lRGL

# =========================
# SOURCES
# =========================
SOURCES := $(shell find $(SRC_DIRS) -name "*.cpp")
OBJECTS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# =========================
# DEFAULT RULE
# =========================
all: $(TARGET) $(GUI_TARGET)

# =========================
# LINK
# =========================
$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# =========================
# COMPILE
# =========================
# -MMD -MP: generate .d dependency files so header changes trigger recompile.
DEPFLAGS = -MMD -MP
-include $(patsubst %.cpp,$(OBJ_DIR)/%.d,$(SOURCES))

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

# =========================
# CREATE DIRS
# =========================
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# =========================
# CLEAN
# =========================
# Removes only compiled artifacts. fit_caches/ and Gamma_fits/ are
# intentionally preserved so accumulated fit results are never lost.
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Full wipe including fit caches and output — use only when you want a
# completely fresh start.
distclean: clean
	rm -rf fit_caches Gamma_fits

# =========================
# RUN
# =========================
run: all
	./$(TARGET)

# Core objects: everything except the CLI main
CORE_SRCS := $(filter-out app/FitGammaAnalyzer.cpp, \
                 $(shell find core fitting tracking models database io analysis_modules -name "*.cpp"))
CORE_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(CORE_SRCS))

# GUI-specific sources
GUI_SRCS     := $(wildcard $(GUI_DIR)/*.cpp)
GUI_OWN_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(GUI_SRCS))

# rootcling dictionary
GUI_DICT_SRC := $(OBJ_DIR)/gui/GammaFitGUIDict.cpp
GUI_DICT_OBJ := $(OBJ_DIR)/gui/GammaFitGUIDict.o

ROOTCLING := $(shell root-config --bindir)/rootcling

$(GUI_DICT_SRC): $(GUI_DIR)/GammaFitGUI.h $(GUI_DIR)/LinkDef.h
	@mkdir -p $(dir $@)
	$(ROOTCLING) -f $@ -D_POSIX_SEM_VALUE_MAX=32767 \
	    "-I$(CURDIR)" "-I$(CURDIR)/app" "-I$(CURDIR)/core" \
	    "-I$(CURDIR)/fitting" "-I$(CURDIR)/tracking" "-I$(CURDIR)/models" \
	    "-I$(CURDIR)/database" "-I$(CURDIR)/io" "-I$(CURDIR)/analysis_modules" \
	    "-I$(CURDIR)/gui" -I$(shell root-config --incdir) $^

$(GUI_DICT_OBJ): $(GUI_DICT_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(patsubst %.cpp,$(OBJ_DIR)/%.d,$(GUI_SRCS))

$(OBJ_DIR)/gui/%.o: $(GUI_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

.PHONY: gui run-gui
gui: $(GUI_TARGET)

run-gui: gui
	ROOT_INCLUDE_PATH=".:./app:./core:./fitting:./tracking:./models:./database:./io:./analysis_modules:./gui" ./$(GUI_TARGET)

$(GUI_TARGET): $(CORE_OBJS) $(GUI_OWN_OBJS) $(GUI_DICT_OBJ) | $(BIN_DIR)
	$(CXX) $^ -o $@ $(GUI_ROOTLIBS)
	cp $(OBJ_DIR)/gui/GammaFitGUIDict_rdict.pcm $(BIN_DIR)/