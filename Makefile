BASEDIR=.
SRCDIR=$(BASEDIR)/src
OBJDIR=$(BASEDIR)/obj

CPP=g++
HOSTFLAGS=-O3 -w -pthread -fopenmp
TARGET=spmm
HEADER=$(SRCDIR)/SpMM.hpp

SRCS=$(notdir $(wildcard $(SRCDIR)/*.cpp))
OBJS=$(SRCS:.cpp=.o) 
OBJECTS = $(patsubst %.o,$(OBJDIR)/%.o,$(OBJS))


all: ${TARGET}

${OBJDIR}/%.o: ${SRCDIR}/%.cpp
	${CPP} ${HOSTFLAGS} -c $< -o $@

${TARGET}: ${OBJECTS}
	${CPP} ${HOSTFLAGS} -o $(TARGET) ${OBJECTS}

# spgemm: ${OBJS}
# 	${CPP} ${HOSTFLAGS} -c -o $@ 

# spgemm: main.o readInput.o compress.o merge.o util.o
# 	$(CPP) $(HOSTFLAGS) -o compress $^

# main.o: main.c $(HEADER)
# 	$(CPP) $(HOSTFLAGS) -c -o $@ $<

# readInput.o: readInput.c $(HEADER)
# 	$(CPP) $(HOSTFLAGS) -c -o $@ $<

# compress.o: compress.c $(HEADER)
# 	$(CPP) $(HOSTFLAGS) -c -o $@ $<

# merge.o: merge.c $(HEADER)
# 	$(CPP) $(HOSTFLAGS) -c -o $@ $<

# util.o: util.c $(HEADER)
# 	$(CPP) $(HOSTFLAGS) -c -o $@ $<


clean:
	rm ${OBJDIR}/*.o ${TARGET}
