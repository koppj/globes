# GLoBES -- General LOng Baseline Experiment Simulator
# (C) 2002 - 2004,  The GLoBES Team
#
# GLoBES is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBES or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


# This is a simple Makefile for the GLoBES examples
#
# Compile example files with ``make example1'' etc.
#
# This Makefile assumes dynamic libraries installed without 
# superuser privileges, where the library is installed to  
# LD_RUN_PATH (set path in your bash configuration file)
# The header file needs to be found in LD_INCLUDE_PATH
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1
#
# For using static libraries, see in the README file
#
# For linking against an installed dynamic library (needs superuser
# privileges), INCLFLAGS becomes obsolete and LDFLAGS simply reads
# -lglobes -lstdc++; the paths do not have to be assigned then

CFLAGS = 
LDFLAGS = ${LD_RUN_PATH}/libglobes.so -lstdc++ 
INCLFLAGS = -I${GLB_INCLUDE_PATH}/ 

BIN = example1 example2 example3 example4
OBJ = example1.o example2.o example3.o example4.o myio.o

example1: example1.o myio.o
	gcc $(CFLAGS) example1.o myio.o -o example1 $(LDFLAGS)

example2: example2.o myio.o
	gcc $(CFLAGS) example2.o myio.o -o example2 $(LDFLAGS)

example3: example3.o myio.o
	gcc $(CFLAGS) example3.o myio.o -o example3 $(LDFLAGS)

example4: example4.o myio.o
	gcc $(CFLAGS) example4.o myio.o -o example4 $(LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) -c $< $(INCLFLAGS)
.PHONY: clean
clean:
	rm $(BIN) $(OBJ)
