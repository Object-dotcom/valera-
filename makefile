CC:= gcc
LIBDIR:= minIni/dev/
CFLAGS:= -I$(LIBDIR)
LDFLAGS:= -lm

NAME:= sitnikov
SRCS = main.c $(LIBDIR)minIni.c
OBJS = main.o $(LIBDIR)minIni.o

.PHONY: all $(NAME) %.o

all: $(NAME)

$(NAME): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

run: $(NAME)
	./$(NAME)

clean:
	rm -f *.o $(LIBDIR)*.o $(NAME)
