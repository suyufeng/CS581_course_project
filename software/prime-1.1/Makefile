### you may change the three options below depending on your favorite compiler ###
# these options mainly for GCC v4.0 or later; previous GCCs cannot specify some options

#### for GCC C++ Compiler
CXX		:= g++# -DMSA_TEST -DMSA_MATRIX
# compiler warning and optimization
# in some cases, `-DCLIMITS' may also be specified
# for GCC v3
CXXFLAGS	:= -Wall -W \
		-Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal \
		-Wpointer-arith -Wswitch-enum -Woverloaded-virtual \
		-O2
# for GCC v4
#CXXFLAGS	:= -Wall -Wextra -Wformat=2 -Wstrict-aliasing=2\
		-Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal \
		-Wpointer-arith -Wswitch-enum -Woverloaded-virtual \
		-O2
# architecture-specific optimization
# If an error has been occurred, use `-march' option instead of `-mtune'
# for Pentium3
#ARCHFLAGS	:= -mtune=pentium3 -mcpu=pentium3
# for Pentium4
#ARCHFLAGS	:= -mtune=pentium4 -mcpu=pentium4
# for G5 Mac
#ARCHFLAGS	:= -mtune=G5 -mcpu=G5

### for Intel C++ Compiler
#CXX		:= icc
# compiler warning and optimization
# in some cases, `-DCLIMITS' may also be specified
#CXXFLAGS	:= -w1 -O2
# architecture-specific optimization
# for Pentium3
#ARCHFLAGS      := -xK
# for Pentium4
#ARCHFLAGS      := -xN
# for Intel64 supporting processors
#ARCHFLAGS      := -xP

# if you would like static binary, comment out the below line
# Note: under Mac OS environment, static binary cannot be created
#STATIC		:= -static


### the following lines should not be modefied ###

LINKER		:= $(CXX) $(STATIC)
LDFLAGS		:=

COMPILE		:= $(CXX) $(CXXFLAGS) $(ARCHFLAGS) -c
LINK		:= $(LINKER) $(LDFLAGS)

STRIP		:= strip

GCC		:= gcc
# Note: You must use GNU sed.
# If not, you have to rewrite make-depend macro.
SED		:= sed

RM		:= rm -f
MV		:= mv -f

SOURCES		:= $(wildcard *.cpp)
OBJS		:= $(subst .cpp,.o,$(SOURCES))

PROGRAM		:= prime

all:
ifneq "$(MAKECMDGOALS)" "clean"
  -include $(subst .cpp,.d,$(SOURCES))
endif

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	@ $(LINK) $^ -o $@
	@ $(STRIP) $@

.PHONY: clean
clean:
	$(RM) *.o *.d $(PROGRAM)

# $(call make-depend,source-file,object-file,depend-file)
define make-depend
  @ $(GCC) -MM -MF $3 -MP -MT $2 $1
endef
# If your gcc does not support the above options,
# you can use the following macro.
# Note that you have to use GNU sed; if not,
# you have to rewrite `sed' part of the macro.
#define make-depend
#  @ $(COMPILE) -M $1 | \
#  $(SED) 's,\($$(notdir $2)\) *:,$$(dir $2) $3: ,' > $3.tmp
#  $(SED) -e 's/#.*//' \
#         -e 's/^[^:]*: *//' \
#         -e 's/ *\\$$$$//'  \
#         -e '/^$$$$/ d'     \
#         -e 's/$$$$/ :/' $3.tmp >> $3.tmp
#  @ $(MV) $3.tmp $3
#endef

%.o: %.cpp
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(COMPILE) $< 

# This target is for those who can not use gcc, GNU make, and/or GNU sed.
.PHONY:	binary
binary: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(ARCHFLAGS) $(STATIC) $(SOURCES) -o $(PROGRAM)
	@ $(STRIP) $(PROGRAM)
