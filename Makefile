# =============================================================================
#  UCMuon Makefile  v2
#  UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
#  MIT License 2026
#
#  Targets:
#    make / make all                  build everything available
#    make local                       OMP-only binaries  (macOS / Linux GUI)
#    make hpc                         MPI+OMP binaries   (Lemaitre4 / CECI)
#    make data-links                  create MUSIC table symlinks in bin/
#    make clean                       remove build/
#    make veryclean                   clean + remove bin/ contents
#    make help                        print this summary
#
#  MUSIC engine (Engine 1) is OPTIONAL.
#  If src/transport/music/music.f is absent the MUSIC targets are skipped.
#  See docs/MUSIC_FILES.md
#
#  Local prerequisites:
#    gfortran with -fopenmp  (brew install gcc  or  apt install gfortran)
#
#  HPC prerequisites (Lemaitre4 / CECI):
#    module load releases/2023b
#    module load foss/2023b
# =============================================================================

FC    = gfortran
MPIFC = mpif90
CC    = gcc
BUILD = build
BIN   = bin

# ---------------------------------------------------------------------------
#  PUMAS backend (optional — requires external/pumas-master/)
# ---------------------------------------------------------------------------
PUMAS_DIR  := external/pumas-master
PUMAS_SRC  := $(PUMAS_DIR)/src/pumas.c
PUMAS_INC  := $(PUMAS_DIR)/include
PUMAS_DRV  := src/transport/pumas/ucmuon_transport_pumas.c
PUMAS_AVAIL := $(shell test -f $(PUMAS_SRC) && echo yes || echo no)

FFLAGS_OMP   = -O2 -fopenmp -Wall -Wno-unused-variable -J$(BUILD) -I$(BUILD)
FFLAGS_PARMA = -O2 -fopenmp -ffree-line-length-none -Wno-implicit-interface \
               -J$(BUILD) -I$(BUILD)
FFLAGS_F77   = -O2 -std=legacy

# ---------------------------------------------------------------------------
#  MUSIC availability check
# ---------------------------------------------------------------------------
MUSIC_SRC := $(wildcard src/transport/music/music.f)
ifeq ($(MUSIC_SRC),)
  MUSIC_AVAIL = no
else
  MUSIC_AVAIL = yes
endif

# ---------------------------------------------------------------------------
#  Source search paths
# ---------------------------------------------------------------------------
vpath %.f90 src/generator src/parma \
            src/transport/music src/transport/bethe_bloch \
            src/common src/converters
vpath %.f   src/common src/transport/music

$(shell mkdir -p $(BUILD) $(BIN))

# ---------------------------------------------------------------------------
#  Object lists
# ---------------------------------------------------------------------------
COMMON_OBJS = \
    $(BUILD)/geom_module_omp.o      \
    $(BUILD)/phits_module_omp.o     \
    $(BUILD)/rng_parallel.o         \
    $(BUILD)/ucmuon_source_module.o \
    $(BUILD)/parma_path_module.o    \
    $(BUILD)/parma_subroutines.o    \
    $(BUILD)/ranlux.o

MUSIC_SUPPORT_OBJS = \
    $(BUILD)/music.o                \
    $(BUILD)/music-crosssections.o  \
    $(BUILD)/ranmar_omp.o           \
    $(BUILD)/rnorml.o               \
    $(BUILD)/corgen.o               \
    $(BUILD)/ranlux_omp.o

GEN_OMP_OBJS   = $(COMMON_OBJS) $(BUILD)/ucmuon_gen_omp.o
MUSIC_OMP_OBJS = $(BUILD)/ucmuon_transport_music_omp.o $(MUSIC_SUPPORT_OBJS)
BB_OMP_OBJS    = $(BUILD)/ucmuon_transport_bb_omp.o $(BUILD)/ranlux_omp.o

UCMUON_GEN_OBJS = $(COMMON_OBJS) $(BUILD)/ucmuon_gen.o
UCMUON_MUS_OBJS = $(BUILD)/ucmuon_transport_music.o $(MUSIC_SUPPORT_OBJS)
UCMUON_BB_OBJS  = $(BUILD)/ucmuon_transport_bb.o $(BUILD)/ranlux_omp.o

# ---------------------------------------------------------------------------
#  Top-level targets
# ---------------------------------------------------------------------------
ifeq ($(MUSIC_AVAIL),yes)
  LOCAL_TARGETS = $(BIN)/ucmuon_gen_omp $(BIN)/ucmuon_transport_music_omp $(BIN)/ucmuon_transport_bb_omp
  HPC_TARGETS   = $(BIN)/ucmuon_gen $(BIN)/ucmuon_transport_music $(BIN)/ucmuon_transport_bb $(BIN)/ucmuon_to_phits
else
  LOCAL_TARGETS = $(BIN)/ucmuon_gen_omp $(BIN)/ucmuon_transport_bb_omp
  HPC_TARGETS   = $(BIN)/ucmuon_gen $(BIN)/ucmuon_transport_bb $(BIN)/ucmuon_to_phits
endif

ifeq ($(PUMAS_AVAIL),yes)
  LOCAL_TARGETS += $(BIN)/ucmuon_transport_pumas
endif

.PHONY: all local hpc data-links clean veryclean help music-status pumas
.PHONY: ucmuon_gen_omp ucmuon_transport_music_omp ucmuon_transport_bb_omp
.PHONY: ucmuon_gen ucmuon_transport_music ucmuon_transport_bb ucmuon_to_phits
.PHONY: ucmuon_transport_pumas

all: local hpc

local: data-links $(LOCAL_TARGETS)
	@$(MAKE) --no-print-directory music-status

hpc: data-links $(HPC_TARGETS)
	@$(MAKE) --no-print-directory music-status

# Convenience aliases (so `make ucmuon_gen_omp` still works)
ucmuon_gen_omp:            $(BIN)/ucmuon_gen_omp
ucmuon_transport_music_omp: $(BIN)/ucmuon_transport_music_omp
ucmuon_transport_bb_omp:   $(BIN)/ucmuon_transport_bb_omp
ucmuon_gen:                $(BIN)/ucmuon_gen
ucmuon_transport_music:    $(BIN)/ucmuon_transport_music
ucmuon_transport_bb:       $(BIN)/ucmuon_transport_bb
ucmuon_to_phits:           $(BIN)/ucmuon_to_phits
ucmuon_transport_pumas:    $(BIN)/ucmuon_transport_pumas
pumas:                     $(BIN)/ucmuon_transport_pumas

music-status:
ifeq ($(MUSIC_AVAIL),no)
	@echo ""
	@echo "  NOTE: Engine 1 (MUSIC) was not built — music.f not found."
	@echo "        Engines 2–6 are fully functional."
	@echo "        See docs/MUSIC_FILES.md to enable Engine 1."
	@echo ""
endif
ifeq ($(PUMAS_AVAIL),no)
	@echo "  NOTE: PUMAS engine not built — external/pumas-master/src/pumas.c not found."
	@echo "        Clone pumas into external/pumas-master/ (see docs/MUSIC_FILES.md)."
	@echo ""
endif

# ---------------------------------------------------------------------------
#  MUSIC data symlinks — created in bin/ (CWD when binaries run)
#  Falls back to cp if ln -sf is unavailable.
# ---------------------------------------------------------------------------
data-links: $(BIN)/music-eloss-rock.dat $(BIN)/music-double-diff-rock.dat \
            $(BIN)/music-cross-sections-rock.dat

$(BIN)/music-eloss-rock.dat: data/music-eloss-rock.dat
	@if [ ! -e $@ ]; then \
	    ln -sf $(abspath $<) $@ 2>/dev/null || cp $< $@; \
	    echo "  LINK  $@ -> $<"; \
	fi

$(BIN)/music-double-diff-rock.dat: data/music-double-diff-rock.dat
	@if [ ! -e $@ ]; then \
	    ln -sf $(abspath $<) $@ 2>/dev/null || cp $< $@; \
	    echo "  LINK  $@ -> $<"; \
	fi

$(BIN)/music-cross-sections-rock.dat: data/music-cross-sections-rock.dat
	@if [ ! -e $@ ]; then \
	    ln -sf $(abspath $<) $@ 2>/dev/null || cp $< $@; \
	    echo "  LINK  $@ -> $<"; \
	fi

# ---------------------------------------------------------------------------
#  Binaries — OMP only (local / GUI)  →  go to bin/
# ---------------------------------------------------------------------------
$(BIN)/ucmuon_gen_omp: $(GEN_OMP_OBJS)
	$(FC) -O2 -fopenmp -o $@ $^ -lm
	@echo "  OK   bin/ucmuon_gen_omp  (OMP, local)"

$(BIN)/ucmuon_transport_bb_omp: $(BB_OMP_OBJS)
	$(FC) -O2 -fopenmp -o $@ $^ -lm
	@echo "  OK   bin/ucmuon_transport_bb_omp  (OMP, local)"

ifeq ($(MUSIC_AVAIL),yes)
$(BIN)/ucmuon_transport_music_omp: data-links $(MUSIC_OMP_OBJS)
	$(FC) -O2 -fopenmp -o $@ $(MUSIC_OMP_OBJS) -lm
	@echo "  OK   bin/ucmuon_transport_music_omp  (OMP, local)"
endif

ifeq ($(PUMAS_AVAIL),yes)
$(BIN)/ucmuon_transport_pumas: $(PUMAS_DRV) $(PUMAS_SRC)
	$(CC) -O2 -o $@ $(PUMAS_DRV) $(PUMAS_SRC) -I$(PUMAS_INC) -lm
	@echo "  OK   bin/ucmuon_transport_pumas  (PUMAS backward/forward MC)"
endif

# ---------------------------------------------------------------------------
#  Binaries — MPI+OMP (HPC)  →  go to bin/
# ---------------------------------------------------------------------------
$(BIN)/ucmuon_gen: $(UCMUON_GEN_OBJS)
	$(MPIFC) -O2 -fopenmp -o $@ $^ -lm
	@echo "  OK   bin/ucmuon_gen  (MPI+OMP)"

$(BIN)/ucmuon_transport_bb: $(UCMUON_BB_OBJS)
	$(MPIFC) -O2 -fopenmp -o $@ $^ -lm
	@echo "  OK   bin/ucmuon_transport_bb  (MPI+OMP)"

$(BIN)/ucmuon_to_phits: $(BUILD)/ucmuon_to_phits.o
	$(FC) -O2 -o $@ $^ -lm
	@echo "  OK   bin/ucmuon_to_phits"

ifeq ($(MUSIC_AVAIL),yes)
$(BIN)/ucmuon_transport_music: data-links $(UCMUON_MUS_OBJS)
	$(MPIFC) -O2 -fopenmp -o $@ $(UCMUON_MUS_OBJS) -lm
	@echo "  OK   bin/ucmuon_transport_music  (MPI+OMP)"
endif

# ---------------------------------------------------------------------------
#  Object compilation (unchanged)
# ---------------------------------------------------------------------------
$(BUILD)/geom_module_omp.o: geom_module.f90
	$(FC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/phits_module_omp.o: phits_module.f90
	$(FC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/rng_parallel.o: rng_parallel.f90
	$(FC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/ucmuon_source_module.o: ucmuon_source_module.f90 $(BUILD)/rng_parallel.o
	$(FC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/parma_path_module.o: parma_path_module.f90
	$(FC) $(FFLAGS_PARMA) -c $< -o $@

$(BUILD)/parma_subroutines.o: parma_subroutines.f90 $(BUILD)/parma_path_module.o
	$(FC) $(FFLAGS_PARMA) -c $< -o $@

$(BUILD)/ucmuon_gen_omp.o: ucmuon_gen_omp.f90 \
    $(BUILD)/geom_module_omp.o $(BUILD)/phits_module_omp.o \
    $(BUILD)/rng_parallel.o    $(BUILD)/ucmuon_source_module.o \
    $(BUILD)/parma_path_module.o
	$(FC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/ucmuon_transport_bb_omp.o: ucmuon_transport_bb_omp.f90
	$(FC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/ucmuon_gen.o: ucmuon_gen.f90 \
    $(BUILD)/geom_module_omp.o $(BUILD)/phits_module_omp.o \
    $(BUILD)/rng_parallel.o    $(BUILD)/ucmuon_source_module.o \
    $(BUILD)/parma_path_module.o
	$(MPIFC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/ucmuon_transport_bb.o: ucmuon_transport_bb.f90
	$(MPIFC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/ucmuon_to_phits.o: ucmuon_to_phits.f90
	$(FC) -O2 -c $< -o $@

$(BUILD)/ranlux.o: ranlux.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD)/ranlux_omp.o: ranlux_omp.f
	$(FC) $(FFLAGS_F77) -fopenmp -c $< -o $@

$(BUILD)/ranmar_omp.o: ranmar_omp.f
	$(FC) $(FFLAGS_F77) -fopenmp -c $< -o $@

$(BUILD)/rnorml.o: rnorml.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD)/corset.o: corset.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD)/corgen.o: corgen.f90
	$(FC) $(FFLAGS_F77) -c $< -o $@

ifeq ($(MUSIC_AVAIL),yes)
$(BUILD)/ucmuon_transport_music_omp.o: ucmuon_transport_music_omp.f90
	$(FC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/ucmuon_transport_music.o: ucmuon_transport_music.f90
	$(MPIFC) $(FFLAGS_OMP) -c $< -o $@

$(BUILD)/music.o: music.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD)/music-crosssections.o: music-crosssections.f
	$(FC) $(FFLAGS_F77) -c $< -o $@
endif

# ---------------------------------------------------------------------------
#  Housekeeping
# ---------------------------------------------------------------------------
help:
	@echo ""
	@echo "  UCMuon v2 — build targets"
	@echo "    make / make all    build everything available"
	@echo "    make local         OMP only  (macOS / GUI)"
	@echo "    make hpc           MPI+OMP   (Lemaitre4 / CECI)"
	@echo "    make data-links    MUSIC table symlinks in bin/"
	@echo "    make clean         remove build/"
	@echo "    make veryclean     clean + remove bin/ contents"
	@echo "  MUSIC engine status: $(MUSIC_AVAIL)"
	@echo ""

clean:
	rm -rf $(BUILD)
	@echo "  build/ removed."

veryclean: clean
	rm -f $(BIN)/ucmuon_gen_omp $(BIN)/ucmuon_transport_music_omp \
	      $(BIN)/ucmuon_transport_bb_omp $(BIN)/ucmuon_gen \
	      $(BIN)/ucmuon_transport_music $(BIN)/ucmuon_transport_bb \
	      $(BIN)/ucmuon_to_phits $(BIN)/ucmuon_transport_pumas \
	      $(BIN)/pumas_*.pumas \
	      $(BIN)/music-eloss-rock.dat $(BIN)/music-double-diff-rock.dat \
	      $(BIN)/music-cross-sections-rock.dat
	@echo "  bin/ contents removed."
