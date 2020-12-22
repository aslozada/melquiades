PROG =	melquiades.x

SRCS =	m_averages.f90 m_boxtype.f90 m_cells.f90 m_configuration.f90 \
	m_constants.f90 m_correlation.f90 m_cutoffs.f90 m_dealloc.f90 \
	m_error.f90 m_fluctuations.f90 m_head.f90 m_init.f90 m_inquire.f90 \
	m_interaction.f90 m_isobaric.f90 m_kind.f90 m_ljones.f90 m_longs.f90 \
	m_metropolis.f90 m_pbc.f90 m_plots.f90 m_precells.f90 m_random.f90 \
	m_read.f90 m_references.f90 m_rotation.f90 m_sampler.f90 m_setup.f90 \
	m_shift.f90 m_simtype.f90 m_store.f90 m_unit.f90 m_uplist.f90 \
	m_zeros.f90 melquiades.f90

OBJS =	m_averages.o m_boxtype.o m_cells.o m_configuration.o m_constants.o \
	m_correlation.o m_cutoffs.o m_dealloc.o m_error.o m_fluctuations.o \
	m_head.o m_init.o m_inquire.o m_interaction.o m_isobaric.o m_kind.o \
	m_ljones.o m_longs.o m_metropolis.o m_pbc.o m_plots.o m_precells.o \
	m_random.o m_read.o m_references.o m_rotation.o m_sampler.o m_setup.o \
	m_shift.o m_simtype.o m_store.o m_unit.o m_uplist.o m_zeros.o \
	melquiades.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = gfortran
F90FLAGS = -Wall -g -O3 -ffpe-trap=zero -fbacktrace -fcheck=all -fbounds-check
LDFLAGS = -Wall -g -O3 -ffpe-trap=zero -fbacktrace -fcheck=all -fbounds-check

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

m_averages.o: m_boxtype.o m_configuration.o m_constants.o m_fluctuations.o \
	m_kind.o m_simtype.o m_unit.o m_zeros.o
m_boxtype.o: m_kind.o
m_cells.o: m_kind.o m_simtype.o
m_configuration.o: m_boxtype.o m_interaction.o m_kind.o m_precells.o \
	m_simtype.o m_zeros.o
m_constants.o: m_kind.o
m_correlation.o: m_kind.o m_simtype.o m_unit.o
m_cutoffs.o: m_boxtype.o m_kind.o m_simtype.o
m_dealloc.o: m_boxtype.o
m_fluctuations.o: m_boxtype.o m_constants.o m_kind.o m_simtype.o m_unit.o \
	m_zeros.o
m_head.o: m_unit.o
m_init.o: m_boxtype.o m_cells.o m_constants.o m_cutoffs.o m_kind.o m_longs.o \
	m_precells.o m_simtype.o m_unit.o m_zeros.o
m_inquire.o: m_unit.o
m_interaction.o: m_boxtype.o m_constants.o m_kind.o m_ljones.o m_pbc.o \
	m_precells.o m_simtype.o m_zeros.o
m_isobaric.o: m_boxtype.o m_configuration.o m_constants.o m_kind.o m_random.o \
	m_simtype.o m_zeros.o
m_ljones.o: m_boxtype.o m_constants.o m_kind.o m_longs.o m_simtype.o \
	m_zeros.o
m_longs.o: m_boxtype.o m_constants.o m_kind.o m_simtype.o m_zeros.o
m_metropolis.o: m_averages.o m_boxtype.o m_configuration.o m_correlation.o \
	m_error.o m_interaction.o m_isobaric.o m_kind.o m_pbc.o m_plots.o \
	m_random.o m_rotation.o m_sampler.o m_shift.o m_simtype.o m_store.o \
	m_unit.o m_uplist.o m_zeros.o
m_pbc.o: m_boxtype.o m_kind.o m_simtype.o
m_plots.o: m_kind.o m_simtype.o m_unit.o
m_precells.o: m_boxtype.o m_constants.o m_cutoffs.o m_kind.o m_simtype.o
m_random.o: m_kind.o
m_read.o: m_boxtype.o m_constants.o m_error.o m_head.o m_inquire.o m_kind.o \
	m_simtype.o m_unit.o
m_rotation.o: m_boxtype.o m_constants.o m_kind.o m_random.o m_simtype.o
m_sampler.o: m_boxtype.o m_configuration.o m_constants.o m_kind.o m_simtype.o \
	m_zeros.o
m_setup.o: m_boxtype.o m_configuration.o m_dealloc.o m_error.o m_init.o \
	m_kind.o m_metropolis.o m_read.o m_simtype.o m_unit.o m_zeros.o
m_shift.o: m_boxtype.o m_kind.o m_simtype.o m_unit.o
m_simtype.o: m_kind.o
m_store.o: m_boxtype.o m_kind.o m_simtype.o m_unit.o
m_uplist.o: m_boxtype.o m_kind.o m_precells.o m_simtype.o
m_zeros.o: m_boxtype.o m_kind.o m_simtype.o
melquiades.o: m_boxtype.o m_configuration.o m_dealloc.o m_error.o m_head.o \
	m_init.o m_kind.o m_metropolis.o m_read.o m_references.o m_setup.o \
	m_simtype.o m_unit.o m_zeros.o
