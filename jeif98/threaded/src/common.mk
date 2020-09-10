
         PROGRAMS = $(BIN)/opflow $(BIN)/uwogamma

all::	$(PROGRAMS)

$(BIN)/opflow:  $(OBD)/opflow.o $(OBD)/tens.o $(OBD)/imio.o $(OBD)/inpo.o \
                $(OBD)/opfl.o $(OBD)/velo.o
		$(RM) $@
		$(CC) $(CFLAGS) $(LIBRARYDIR) -o $@ $(OBD)/opflow.o \
                $(OBD)/tens.o $(OBD)/imio.o $(OBD)/inpo.o $(OBD)/opfl.o \
                $(OBD)/velo.o -lm
$(BIN)/uwogamma: $(OBD)/uwogamma.o $(OBD)/tens.o $(OBD)/imio.o
		$(RM) $@
		$(CC) $(CFLAGS) $(LIBRARYDIR) -o $@ $(OBD)/uwogamma.o \
                $(OBD)/tens.o $(OBD)/imio.o -lm

$(OBD)/%.o: %.c
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c $<

$(OBD)/opflow.o: opflow.c glob.h tens.h imio.h inpo.h opfl.h velo.h
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c opflow.c
$(OBD)/uwogamma.o: uwogamma.c glob.h tens.h imio.h
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c uwogamma.c

$(OBD)/tens.o:  tens.c glob.h tens.h
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c tens.c
$(OBD)/imio.o:  imio.c glob.h tens.h imio.h
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c imio.c
$(OBD)/inpo.o:  inpo.c glob.h tens.h imio.h inpo.h velo.h
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c inpo.c
$(OBD)/velo.o:  velo.c glob.h tens.h imio.h inpo.h velo.h
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c velo.c
$(OBD)/opfl.o:  opfl.c glob.h tens.h imio.h inpo.h opfl.h
		$(RM) $@
		$(CC) $(CFLAGS) -o $@ -c opfl.c
