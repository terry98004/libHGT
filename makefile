CC = gcc
CFLAGS = -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c
AR = ar
ARFLAGS = rcs
SRCS = hgtInit.c ThetaOfT.c GramAtN.c GramNearT.c RSbuildcoeff.c RSremainder.c RSmainTerm.c HardyZcalc.c 
OBJS = $(SRCS:.c=.o)
DEPS = hgt.h
TARGET = libhgt.a

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
