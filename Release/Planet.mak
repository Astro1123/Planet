#-------------------------------------------------------
# BCC Developer 1.2.21
# Copyright (C) 2003 jun_miura@hi-ho.ne.jp
#-------------------------------------------------------
.autodepend
CC=bcc32
RC=brc32
CFLAG=-W  -3 -O2 -w- -AT -pc -H- -k -b  
OUTDIR=-nRelease
CINCS=
TARGET=Release\Planet.exe
SRC1=C:\Users\Adans\Desktop\Planet\Planet.cpp
OBJ1=Release\Planet.obj
SRC2=C:\Users\Adans\Desktop\Planet\calc.cpp
OBJ2=Release\calc.obj

TARGET: $(TARGET)

$(TARGET): $(OBJ1) $(OBJ2)
    $(CC) $(CFLAG) -e$(TARGET) $(OBJ1) $(OBJ2)

$(OBJ1): $(SRC1)
    $(CC) $(CFLAG) $(OUTDIR) $(CINCS) -c $(SRC1)

$(OBJ2): $(SRC2)
    $(CC) $(CFLAG) $(OUTDIR) $(CINCS) -c $(SRC2)
