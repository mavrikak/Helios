##### =================================================================================================
##### Toolchains
##### =================================================================================================
# Conda sets $FC and $CXX automatically. 
# ?= ensures we use Conda's compilers if active, or system defaults if not.
FC  ?= gfortran
CXX ?= g++

# Point custom variables to the resolved CXX compiler
CPP = $(CXX) -std=c++11
LL  = $(CXX)

##### =================================================================================================
##### Warnings & Optimization
##### =================================================================================================
WARN_FLAGS        = -W -Wall -Wno-lto-type-mismatch -Wno-sign-compare
# Default: avoid cross-language LTO surprises (enable with `make lto` if desired)
OPT_FLAGS         = -O3 -m64 -fno-strict-aliasing
OPT_FLAGS_AMOS    = -O3 -m64 -fno-lto -fno-strict-aliasing
OPT_FLAGS_QUAD    = -O3 -m64 -fno-lto -fno-strict-aliasing
OPT_FLAGS_REGRID  = -O3 -m64 -fno-lto -fno-strict-aliasing
STRIP_FLAGS       = -s

# Pick up the active conda env path
CONDA_PREFIX ?= $(shell echo $$CONDA_PREFIX)

# Keep just one LIB_FLAGS line, augmented with conda paths:
LIB_FLAGS = -lopenblas -lgfortran -lpthread
ifdef CONDA_PREFIX
LIB_FLAGS += -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
endif

##### =================================================================================================
##### Layout
##### =================================================================================================
BUILD   := build
APP_DIR := apps
OBJ_DIR := $(BUILD)/obj
LIB_DIR := $(BUILD)/lib
MOD_DIR := $(BUILD)/mod

SRC_DIR := src
INC_DIR := include
APP_DIR := apps
TP_DIR  := third_party

# Fortran modules: write with -J, search with -I
FCMOD    := -J$(MOD_DIR)
FCSEARCH := -I$(MOD_DIR)

# Ensure dirs exist
$(shell mkdir -p $(APP_DIR) $(OBJ_DIR) $(LIB_DIR) $(MOD_DIR))

##### =================================================================================================
##### Sources
##### =================================================================================================
# Apps (mains)
APPS_SRC := $(APP_DIR)/main.cpp \
            $(APP_DIR)/mainpp.cpp \
            $(APP_DIR)/casimir.cpp \
            $(APP_DIR)/maintable.cpp

# Core C++
CPP_SRC  := $(wildcard $(SRC_DIR)/*.cpp)

# Fortran third-party
REGRID_F90 := $(TP_DIR)/regridpack/regridpack_module.f90 \
              $(TP_DIR)/regridpack/regridpack_c.f90
QUAD_F     := $(TP_DIR)/quadpack/dqage.f
LEGACY_F   := $(TP_DIR)/legacy_f/cg.f
AMOS_DIR   := $(TP_DIR)/amos
AMOS_F     := $(wildcard $(AMOS_DIR)/zbesj.f) \
              $(wildcard $(AMOS_DIR)/zbesh.f) \
              $(wildcard $(AMOS_DIR)/zbesy.f) \
              $(wildcard $(AMOS_DIR)/d1mach.f) \
              $(wildcard $(AMOS_DIR)/i1mach.f)

##### =================================================================================================
##### Derived names
##### =================================================================================================
CPP_OBJ     := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_SRC))
APP_OBJ     := $(patsubst $(APP_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(APPS_SRC))

REGRID_OBJ  := $(OBJ_DIR)/regridpack_module.o \
               $(OBJ_DIR)/regridpack_c.o
QUAD_OBJ    := $(patsubst $(TP_DIR)/quadpack/%.f,$(OBJ_DIR)/quadpack_%.o,$(QUAD_F))
LEGACY_OBJ  := $(patsubst $(TP_DIR)/legacy_f/%.f,$(OBJ_DIR)/legacy_%.o,$(LEGACY_F))
AMOS_OBJ    := $(patsubst $(AMOS_DIR)/%.f,$(OBJ_DIR)/amos_%.o,$(AMOS_F))


LIBSIM      := $(LIB_DIR)/libsimcore.a
LIBREGRID   := $(LIB_DIR)/libregridpack.a
LIBQUAD     := $(LIB_DIR)/libquadpack.a
LIBLEGACY   := $(LIB_DIR)/liblegacy.a

SIENANO     := $(APP_DIR)/SIENano
SIENANOPP   := $(APP_DIR)/SIENanoPP
CASIMIR     := $(APP_DIR)/Casimir
SIETABLE    := $(APP_DIR)/SIETable

##### =================================================================================================
##### Flags
##### =================================================================================================
CCFLAGS        = $(OPT_FLAGS) -MMD -MP -I$(INC_DIR) -I$(SRC_DIR) # $(WARN_FLAGS)
FCFLAGS        = $(OPT_FLAGS) -fimplicit-none
FCFLAGS_AMOS   = $(OPT_FLAGS_AMOS) -std=legacy -ffixed-form -Wno-maybe-uninitialized -Wno-uninitialized
FCFLAGS_QUAD   = $(OPT_FLAGS_QUAD) -std=legacy -ffixed-form -Wno-maybe-uninitialized -Wno-uninitialized
FCFLAGS_REGRID = $(OPT_FLAGS_REGRID) -fimplicit-none

LDFLAGS = $(STRIP_FLAGS) -Wl,--no-as-needed $(LIB_FLAGS) # $(WARN_FLAGS)

##### =================================================================================================
##### Phony
##### =================================================================================================
.PHONY: all clean debug gdb lto decolourised

all: $(SIENANO) $(SIENANOPP) $(CASIMIR)

table: $(SIETABLE)

##### =================================================================================================
##### Link executables
##### =================================================================================================
$(SIENANO): $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(OBJ_DIR)/main.o
	@echo "Link  $@"
	@$(LL) -o $@ $(OBJ_DIR)/main.o $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(LDFLAGS)

$(SIENANOPP): $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(OBJ_DIR)/mainpp.o
	@echo "Link  $@"
	@$(LL) -o $@ $(OBJ_DIR)/mainpp.o $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(LDFLAGS)

$(CASIMIR): $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(OBJ_DIR)/casimir.o
	@echo "Link  $@"
	@$(LL) -o $@ $(OBJ_DIR)/casimir.o $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(LDFLAGS)

$(SIETABLE): $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(OBJ_DIR)/maintable.o
	@echo "Link  $@"
	@$(LL) -o $@ $(OBJ_DIR)/maintable.o $(LIBSIM) $(LIBREGRID) $(LIBQUAD) $(LIBLEGACY) $(LDFLAGS)

##### =================================================================================================
##### Static libraries
##### =================================================================================================
$(LIBSIM): $(CPP_OBJ) | $(LIB_DIR)
	@echo "Ar    $@"
	@ar rcs $@ $(CPP_OBJ)

$(LIBREGRID): $(REGRID_OBJ) | $(LIB_DIR)
	@echo "Ar    $@"
	@ar rcs $@ $(REGRID_OBJ)

$(LIBQUAD): $(QUAD_OBJ) | $(LIB_DIR)
	@echo "Ar    $@"
	@ar rcs $@ $(QUAD_OBJ)

$(LIBLEGACY): $(LEGACY_OBJ) $(AMOS_OBJ) | $(LIB_DIR)
	@echo "Ar    $@"
	@ar rcs $@ $(LEGACY_OBJ) $(AMOS_OBJ)

##### =================================================================================================
##### Compile rules — C++
##### =================================================================================================
$(OBJ_DIR)/%.o: $(APP_DIR)/%.cpp
	@echo "C++   $<"
	@$(CPP) -c $< -o $@ $(CCFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "C++   $<"
	@$(CPP) -c $< -o $@ $(CCFLAGS)

##### =================================================================================================
##### Compile rules — Fortran
##### =================================================================================================
# Build module first (writes .mod to MOD_DIR), then the user wrapper.
$(OBJ_DIR)/regridpack_module.o: $(TP_DIR)/regridpack/regridpack_module.f90 | $(MOD_DIR)
	@echo "F90   $<"
	@$(FC) -c $< -o $@ $(FCFLAGS_REGRID) $(FCMOD)

$(OBJ_DIR)/regridpack_c.o: $(TP_DIR)/regridpack/regridpack_c.f90 $(OBJ_DIR)/regridpack_module.o | $(MOD_DIR)
	@echo "F90   $<"
	@$(FC) -c $< -o $@ $(FCFLAGS_REGRID) $(FCMOD) $(FCSEARCH)

# QUADPACK (dqage)
$(OBJ_DIR)/quadpack_%.o: $(TP_DIR)/quadpack/%.f
	@echo "F77   $<"
	@$(FC) -c $< -o $@ $(FCFLAGS_QUAD) $(FCMOD) $(FCSEARCH)

# Legacy fixed-form (cg)
$(OBJ_DIR)/legacy_%.o: $(TP_DIR)/legacy_f/%.f
	@echo "F77   $<"
	@$(FC) -c $< -o $@ $(FCFLAGS) $(FCMOD) $(FCSEARCH)

# AMOS (zbesj/zbesy/zbesh/d1mach/i1mach)
$(OBJ_DIR)/amos_%.o: $(TP_DIR)/amos/%.f
	@echo "F77   $<"
	@$(FC) -c $< -o $@ $(FCFLAGS_AMOS) $(FCMOD) $(FCSEARCH)

##### =================================================================================================
##### Convenience modes
##### =================================================================================================
profile: OPT_FLAGS = -pg -O3 -m64 -fno-strict-aliasing
profile: STRIP_FLAGS =
profile: all

debug: CCFLAGS += -DDEBUG -g
debug: STRIP_FLAGS =
debug: all

gdb: OPT_FLAGS = -g
gdb: STRIP_FLAGS =
gdb: all

# Optional: try C++ LTO (Fortran stays -fno-lto for stability)
lto: CCFLAGS += -flto
lto: all

##### =================================================================================================
##### Clean
##### =================================================================================================
clean:
	@echo "Clean"
	@rm -rf $(BUILD)

##### =================================================================================================
##### Auto-deps for C++
##### =================================================================================================
-include $(CPP_OBJ:.o=.d) $(APP_OBJ:.o=.d)
