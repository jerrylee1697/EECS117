.DEFAULT_GOAL := all

NVCC = nvcc

TARGETS = naive$(EXEEXT) stride$(EXEEXT) sequential$(EXEEXT) first_add$(EXEEXT) unroll$(EXEEXT) multiple$(EXEEXT) transpose$(EXEEXT)

all: $(TARGETS)

SRCS_COMMON = timer.c

DISTFILES += $(SRCS_COMMON) $(DEPS_COMMON)

naive$(EXEEXT): stride.cu $(SRCS_COMMON)
	$(NVCC) -o $@ naive.cu $(SRCS_COMMON) 

stride$(EXEEXT): stride.cu $(SRCS_COMMON)
	$(NVCC) -o $@ stride.cu $(SRCS_COMMON)

sequential$(EXEEXT): sequential.cu $(SRCS_COMMON)
	$(NVCC) -o $@ sequential.cu $(SRCS_COMMON)

first_add$(EXEEXT): first_add.cu $(SRCS_COMMON)
	$(NVCC) -o $@ first_add.cu $(SRCS_COMMON)

unroll$(EXEEXT): unroll.cu $(SRCS_COMMON)
	$(NVCC) -o $@ unroll.cu $(SRCS_COMMON)

multiple$(EXEEXT): multiple.cu $(SRCS_COMMON)
	$(NVCC) -o $@ multiple.cu $(SRCS_COMMON)

transpose$(EXEEXT): transpose.cu $(SRS_COMMON)
	$(NVCC) -o $@ transpose.cu $(SRCS_COMMON)
clean:
	rm -f $(TARGETS) 

# eof
